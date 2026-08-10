"""Convert the ultra-wide GSE255460 count table to partitioned disk CSR.

The input is a genes-by-cells dense text member inside a tar.gz archive.  The
output remains genes-by-cells but stores only non-zero integer values, split by
the technical ``ID`` partitions in the metadata.  Each partition contains
little-endian int32 CSR arrays plus explicit feature and barcode tables.

The converter validates every source value indirectly through exact, all-cell
nCount/nFeature agreement with the published metadata.  It writes to a
temporary sibling directory and publishes ``manifest.json`` only after all
validation succeeds.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np


FORMAT_NAME = "gse255460_partitioned_csr"
FORMAT_VERSION = 1
SCRIPT_VERSION = "1.0.0"
INT32_MAX = np.iinfo(np.int32).max


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert GSE255460 sc_counts.txt to validated partitioned disk CSR."
        )
    )
    parser.add_argument("--archive", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--member", default="sc_counts.txt")
    parser.add_argument("--batch-column", default="ID")
    parser.add_argument("--cell-column", default="X")
    parser.add_argument("--count-column", default="nCount_RNA")
    parser.add_argument("--feature-column", default="nFeature_RNA")
    parser.add_argument("--buffer-rows", type=int, default=64)
    parser.add_argument("--progress-every", type=int, default=1000)
    return parser.parse_args()


def file_signature(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.resolve()).replace("\\", "/"),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def clean_partition(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    if not cleaned:
        raise ValueError("batch ID cannot be converted to a safe directory name")
    return cleaned


def make_r_name(value: str) -> str:
    """Implement the make.names behavior needed by these ASCII cell IDs."""
    cleaned = re.sub(r"[^A-Za-z0-9._]", ".", value)
    if not re.match(r"^[A-Za-z]", cleaned) and not re.match(
        r"^\.[^0-9]", cleaned
    ):
        cleaned = "X" + cleaned
    if cleaned in {
        "if",
        "else",
        "repeat",
        "while",
        "function",
        "for",
        "in",
        "next",
        "break",
        "TRUE",
        "FALSE",
        "NULL",
        "Inf",
        "NaN",
        "NA",
        "NA_integer_",
        "NA_real_",
        "NA_complex_",
        "NA_character_",
    }:
        cleaned += "."
    return cleaned


def parse_nonnegative_integer(value: str, column: str, row: int) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise ValueError(
            f"metadata row {row} has non-integer {column}: {value!r}"
        ) from error
    if parsed < 0:
        raise ValueError(f"metadata row {row} has negative {column}")
    return parsed


def read_metadata(
    path: Path,
    cell_column: str,
    batch_column: str,
    count_column: str,
    feature_column: str,
) -> dict[str, object]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    cells: list[str] = []
    batches: list[str] = []
    expected_counts: list[int] = []
    expected_features: list[int] = []
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            cell_column,
            batch_column,
            count_column,
            feature_column,
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "metadata is missing required columns: "
                + ", ".join(sorted(missing))
            )
        for row_number, row in enumerate(reader, start=2):
            cell = (row[cell_column] or "").strip()
            batch = (row[batch_column] or "").strip()
            if not cell or not batch:
                raise ValueError(
                    f"metadata row {row_number} has an empty cell or batch ID"
                )
            cells.append(cell)
            batches.append(batch)
            expected_counts.append(
                parse_nonnegative_integer(
                    row[count_column], count_column, row_number
                )
            )
            expected_features.append(
                parse_nonnegative_integer(
                    row[feature_column], feature_column, row_number
                )
            )
    if not cells:
        raise ValueError("metadata has no cells")
    if len(set(cells)) != len(cells):
        raise ValueError("metadata cell IDs are duplicated")

    batch_order = list(dict.fromkeys(batches))
    positions = {
        batch: np.asarray(
            [index for index, value in enumerate(batches) if value == batch],
            dtype=np.int32,
        )
        for batch in batch_order
    }
    safe_names = [clean_partition(batch) for batch in batch_order]
    if len(set(safe_names)) != len(safe_names):
        raise ValueError("batch IDs collide after filesystem sanitization")

    return {
        "cells": cells,
        "matrix_cells": [make_r_name(cell) for cell in cells],
        "batches": batches,
        "batch_order": batch_order,
        "batch_positions": positions,
        "safe_names": dict(zip(batch_order, safe_names)),
        "expected_counts": np.asarray(expected_counts, dtype=np.int64),
        "expected_features": np.asarray(expected_features, dtype=np.int64),
    }


class HashedBinaryWriter:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.handle = path.open("wb")
        self.digest = hashlib.sha256()
        self.size = 0

    def write_array(self, values: np.ndarray) -> None:
        payload = values.tobytes(order="C")
        self.handle.write(payload)
        self.digest.update(payload)
        self.size += len(payload)

    def close(self) -> None:
        self.handle.close()

    @property
    def sha256(self) -> str:
        return self.digest.hexdigest()


def existing_bundle_matches(
    output_dir: Path,
    archive_signature: dict[str, object],
    metadata_signature: dict[str, object],
) -> bool:
    manifest_path = output_dir / "manifest.json"
    if not manifest_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    if (
        manifest.get("format") != FORMAT_NAME
        or manifest.get("format_version") != FORMAT_VERSION
        or manifest.get("script_version") != SCRIPT_VERSION
        or manifest.get("archive_signature") != archive_signature
        or manifest.get("metadata_signature") != metadata_signature
    ):
        return False
    for partition in manifest.get("partitions", []):
        for key in ("data_file", "indices_file", "indptr_file", "barcodes_file"):
            if not (output_dir / partition[key]).is_file():
                return False
    return (output_dir / manifest.get("features_file", "")).is_file()


def write_text_with_sha256(path: Path, lines: list[str]) -> tuple[int, str]:
    payload = "".join(lines).encode("utf-8")
    path.write_bytes(payload)
    return len(payload), hashlib.sha256(payload).hexdigest()


def convert(args: argparse.Namespace) -> None:
    started = time.time()
    archive = Path(args.archive).resolve()
    metadata_path = Path(args.metadata).resolve()
    output_dir = Path(args.output_dir).resolve()
    if not archive.is_file():
        raise FileNotFoundError(f"archive does not exist: {archive}")
    if not metadata_path.is_file():
        raise FileNotFoundError(f"metadata does not exist: {metadata_path}")
    if args.buffer_rows < 1:
        raise ValueError("--buffer-rows must be at least 1")

    archive_signature = file_signature(archive)
    metadata_signature = file_signature(metadata_path)
    if existing_bundle_matches(
        output_dir, archive_signature, metadata_signature
    ):
        print(f"Reusing validated sparse bundle: {output_dir}")
        return

    metadata = read_metadata(
        metadata_path,
        args.cell_column,
        args.batch_column,
        args.count_column,
        args.feature_column,
    )
    cells = metadata["cells"]
    matrix_cells = metadata["matrix_cells"]
    batch_order = metadata["batch_order"]
    batch_positions = metadata["batch_positions"]
    safe_names = metadata["safe_names"]
    n_cells = len(cells)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(
        tempfile.mkdtemp(
            prefix=output_dir.name + ".building-",
            dir=str(output_dir.parent),
        )
    )
    writers: dict[str, dict[str, HashedBinaryWriter]] = {}
    buffers: dict[str, dict[str, list[np.ndarray]]] = {}
    row_nonzero: dict[str, list[int]] = {}
    partition_dirs: dict[str, Path] = {}
    features: list[str] = []
    feature_set: set[str] = set()
    observed_counts = np.zeros(n_cells, dtype=np.int64)
    observed_features = np.zeros(n_cells, dtype=np.int32)
    total_nonzero = 0

    try:
        for batch in batch_order:
            partition_dir = temporary_dir / safe_names[batch]
            partition_dir.mkdir()
            partition_dirs[batch] = partition_dir
            writers[batch] = {
                "data": HashedBinaryWriter(partition_dir / "data.i32"),
                "indices": HashedBinaryWriter(partition_dir / "indices.i32"),
            }
            buffers[batch] = {"data": [], "indices": []}
            row_nonzero[batch] = []

            positions = batch_positions[batch]
            barcode_lines = ["cell_id\tmatrix_cell_id\tglobal_index_1based\n"]
            barcode_lines.extend(
                f"{cells[index]}\t{matrix_cells[index]}\t{index + 1}\n"
                for index in positions
            )
            write_text_with_sha256(
                partition_dir / "barcodes.tsv", barcode_lines
            )

        def flush_buffers() -> None:
            for batch in batch_order:
                if buffers[batch]["data"]:
                    writers[batch]["data"].write_array(
                        np.concatenate(buffers[batch]["data"]).astype(
                            "<i4", copy=False
                        )
                    )
                    writers[batch]["indices"].write_array(
                        np.concatenate(buffers[batch]["indices"]).astype(
                            "<i4", copy=False
                        )
                    )
                    buffers[batch]["data"].clear()
                    buffers[batch]["indices"].clear()

        with tarfile.open(archive, mode="r:*") as tar:
            member = tar.getmember(args.member)
            stream = tar.extractfile(member)
            if stream is None:
                raise RuntimeError("could not open count archive member")
            header = stream.readline().rstrip(b"\r\n")
            header_cells = header.decode("utf-8").split("\t")
            if header_cells != matrix_cells:
                mismatch = next(
                    (
                        index
                        for index, pair in enumerate(
                            zip(header_cells, matrix_cells)
                        )
                        if pair[0] != pair[1]
                    ),
                    min(len(header_cells), len(matrix_cells)),
                )
                raise ValueError(
                    "matrix header does not exactly match make.names(metadata "
                    f"cell IDs); first mismatch index is {mismatch + 1}"
                )

            for row_index, raw_line in enumerate(stream, start=1):
                line = raw_line.rstrip(b"\r\n")
                if not line:
                    raise ValueError(f"empty matrix row at gene {row_index}")
                first_tab = line.find(b"\t")
                if first_tab <= 0:
                    raise ValueError(
                        f"matrix row {row_index} has no gene/value separator"
                    )
                gene = line[:first_tab].decode("utf-8")
                if not gene or gene in feature_set:
                    raise ValueError(
                        f"matrix row {row_index} has an empty/duplicate gene"
                    )
                values = np.fromstring(
                    line[first_tab + 1 :],
                    dtype=np.int64,
                    sep="\t",
                )
                if values.size != n_cells:
                    raise ValueError(
                        f"matrix row {row_index} has {values.size} values; "
                        f"expected {n_cells}"
                    )
                if values.min(initial=0) < 0:
                    raise ValueError(
                        f"matrix row {row_index} contains a negative count"
                    )
                if values.max(initial=0) > INT32_MAX:
                    raise ValueError(
                        f"matrix row {row_index} exceeds int32 count range"
                    )

                features.append(gene)
                feature_set.add(gene)
                observed_counts += values
                observed_features += values > 0

                for batch in batch_order:
                    subset = values[batch_positions[batch]]
                    nonzero_indices = np.flatnonzero(subset).astype(
                        np.int32, copy=False
                    )
                    nonzero_values = subset[nonzero_indices].astype(
                        np.int32, copy=False
                    )
                    row_nonzero[batch].append(int(nonzero_indices.size))
                    total_nonzero += int(nonzero_indices.size)
                    if nonzero_indices.size:
                        buffers[batch]["indices"].append(nonzero_indices)
                        buffers[batch]["data"].append(nonzero_values)

                if row_index % args.buffer_rows == 0:
                    flush_buffers()
                if row_index % args.progress_every == 0:
                    elapsed = time.time() - started
                    print(
                        "GSE255460 sparse conversion: "
                        f"{row_index:,} genes, {total_nonzero:,} nonzeros, "
                        f"{elapsed:.1f} seconds"
                    )
                    sys.stdout.flush()

        flush_buffers()
        for batch in batch_order:
            writers[batch]["data"].close()
            writers[batch]["indices"].close()

        expected_counts = metadata["expected_counts"]
        expected_features = metadata["expected_features"]
        count_mismatch = np.flatnonzero(observed_counts != expected_counts)
        feature_mismatch = np.flatnonzero(
            observed_features.astype(np.int64) != expected_features
        )
        if count_mismatch.size or feature_mismatch.size:
            examples = sorted(
                set(count_mismatch[:5].tolist())
                | set(feature_mismatch[:5].tolist())
            )
            detail = "; ".join(
                (
                    f"{cells[index]} observed counts/features "
                    f"{observed_counts[index]}/{observed_features[index]} "
                    f"expected {expected_counts[index]}/"
                    f"{expected_features[index]}"
                )
                for index in examples
            )
            raise ValueError(
                "all-cell metadata validation failed: " + detail
            )
        if total_nonzero != int(expected_features.sum()):
            raise ValueError(
                "total sparse nonzero count does not equal metadata "
                "nFeature_RNA sum"
            )

        feature_lines = ["gene_index_1based\tgene_id\n"]
        feature_lines.extend(
            f"{index}\t{gene}\n"
            for index, gene in enumerate(features, start=1)
        )
        feature_size, feature_sha = write_text_with_sha256(
            temporary_dir / "features.tsv", feature_lines
        )

        partitions: list[dict[str, object]] = []
        for batch in batch_order:
            row_counts = np.asarray(row_nonzero[batch], dtype=np.int64)
            indptr = np.empty(len(features) + 1, dtype=np.int64)
            indptr[0] = 0
            np.cumsum(row_counts, out=indptr[1:])
            if indptr[-1] > INT32_MAX:
                raise ValueError(
                    f"partition {batch} CSR pointer exceeds int32 range"
                )
            indptr_path = partition_dirs[batch] / "indptr.i32"
            indptr_payload = indptr.astype("<i4", copy=False).tobytes()
            indptr_path.write_bytes(indptr_payload)
            barcode_path = partition_dirs[batch] / "barcodes.tsv"
            positions = batch_positions[batch]
            expected_partition_nonzero = int(
                expected_features[positions].sum()
            )
            if int(indptr[-1]) != expected_partition_nonzero:
                raise ValueError(
                    f"partition {batch} nonzero count disagrees with metadata"
                )
            relative_dir = safe_names[batch]
            partitions.append(
                {
                    "partition_id": batch,
                    "safe_id": relative_dir,
                    "n_cells": int(positions.size),
                    "n_genes": len(features),
                    "nonzero": int(indptr[-1]),
                    "global_index_min_1based": int(positions.min()) + 1,
                    "global_index_max_1based": int(positions.max()) + 1,
                    "global_indices_contiguous": bool(
                        np.all(np.diff(positions) == 1)
                    ),
                    "data_file": f"{relative_dir}/data.i32",
                    "data_bytes": writers[batch]["data"].size,
                    "data_sha256": writers[batch]["data"].sha256,
                    "indices_file": f"{relative_dir}/indices.i32",
                    "indices_bytes": writers[batch]["indices"].size,
                    "indices_sha256": writers[batch]["indices"].sha256,
                    "indptr_file": f"{relative_dir}/indptr.i32",
                    "indptr_bytes": len(indptr_payload),
                    "indptr_sha256": hashlib.sha256(
                        indptr_payload
                    ).hexdigest(),
                    "barcodes_file": f"{relative_dir}/barcodes.tsv",
                    "barcodes_bytes": barcode_path.stat().st_size,
                    "barcodes_sha256": sha256_file(barcode_path),
                }
            )

        manifest = {
            "format": FORMAT_NAME,
            "format_version": FORMAT_VERSION,
            "script_version": SCRIPT_VERSION,
            "created_at": time.strftime(
                "%Y-%m-%dT%H:%M:%S%z", time.localtime()
            ),
            "archive_signature": archive_signature,
            "archive_sha256": sha256_file(archive),
            "metadata_signature": metadata_signature,
            "metadata_sha256": sha256_file(metadata_path),
            "matrix_member": args.member,
            "matrix_member_size": int(member.size),
            "orientation": "genes_by_cells",
            "sparse_encoding": "CSR",
            "index_base": 0,
            "value_dtype": "little_endian_int32",
            "index_dtype": "little_endian_int32",
            "pointer_dtype": "little_endian_int32",
            "n_genes": len(features),
            "n_cells": n_cells,
            "n_partitions": len(partitions),
            "nonzero": int(total_nonzero),
            "total_counts": int(observed_counts.sum()),
            "all_cells_nCount_exact": True,
            "all_cells_nFeature_exact": True,
            "header_mapping": (
                "identical(matrix_header, make.names(metadata_cell_id))"
            ),
            "header_mapping_exact": True,
            "features_unique": len(feature_set) == len(features),
            "features_file": "features.tsv",
            "features_bytes": feature_size,
            "features_sha256": feature_sha,
            "batch_column": args.batch_column,
            "cell_column": args.cell_column,
            "count_column": args.count_column,
            "feature_column": args.feature_column,
            "elapsed_seconds": time.time() - started,
            "partitions": partitions,
        }
        manifest_path = temporary_dir / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        if output_dir.exists():
            stale = output_dir.with_name(
                output_dir.name
                + ".stale-"
                + time.strftime("%Y%m%d-%H%M%S", time.localtime())
            )
            output_dir.rename(stale)
            print(f"Moved stale generated bundle to: {stale}")
        # Some Windows antivirus/indexing configurations transiently deny an
        # atomic directory rename even after every file handle is closed.
        # shutil.move retains the fast rename path and falls back to a
        # copy-and-remove operation within this generated cache directory.
        shutil.move(str(temporary_dir), str(output_dir))
        print(
            "GSE255460 sparse conversion completed: "
            f"{len(features):,} genes, {n_cells:,} cells, "
            f"{total_nonzero:,} nonzeros, "
            f"{time.time() - started:.1f} seconds"
        )
        print(f"Validated sparse manifest: {output_dir / 'manifest.json'}")
    except Exception:
        for batch_writers in writers.values():
            for writer in batch_writers.values():
                if not writer.handle.closed:
                    writer.close()
        print(
            f"Conversion failed; recoverable partial output retained at "
            f"{temporary_dir}",
            file=sys.stderr,
        )
        raise


def main() -> None:
    args = parse_args()
    convert(args)


if __name__ == "__main__":
    main()
