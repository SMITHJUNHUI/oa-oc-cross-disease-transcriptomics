from __future__ import print_function

import argparse
import hashlib
import json
import os
import sys
import tarfile
import time


def parse_args():
    parser = argparse.ArgumentParser(
        description="Stream-validate the ultra-wide GSE255460 count member."
    )
    parser.add_argument("--archive", required=True)
    parser.add_argument("--member", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--spotcheck-cells", type=int, default=5)
    return parser.parse_args()


def parse_first_values(line, count):
    first_tab = line.find(b"\t")
    if first_tab <= 0:
        raise ValueError("row has no gene/value separator")
    gene = line[:first_tab]
    values = []
    start = first_tab + 1
    for _ in range(count):
        end = line.find(b"\t", start)
        if end < 0:
            end = len(line)
        token = line[start:end]
        if not token:
            raise ValueError("empty spot-check count")
        try:
            value = int(token)
        except ValueError:
            raise ValueError("spot-check count is not an integer")
        if value < 0:
            raise ValueError("spot-check count is negative")
        values.append(value)
        start = end + 1
    return gene, values


def main():
    args = parse_args()
    started = time.time()
    spotcheck_cells = max(1, args.spotcheck_cells)
    observed_counts = [0] * spotcheck_cells
    observed_features = [0] * spotcheck_cells
    genes = set()
    row_count = 0
    field_mismatches = 0
    first_mismatch_row = None

    with tarfile.open(args.archive, mode="r:*") as archive:
        member = archive.getmember(args.member)
        stream = archive.extractfile(member)
        if stream is None:
            raise RuntimeError("could not open archive member")
        header = stream.readline().rstrip(b"\r\n")
        if not header:
            raise RuntimeError("count member has no header")
        n_cells = header.count(b"\t") + 1
        expected_fields = n_cells + 1
        header_md5 = hashlib.md5(header).hexdigest()

        for raw_line in stream:
            line = raw_line.rstrip(b"\r\n")
            if not line:
                raise RuntimeError(
                    "empty matrix row at gene row {0}".format(row_count + 1)
                )
            row_count += 1
            observed_fields = line.count(b"\t") + 1
            if observed_fields != expected_fields:
                field_mismatches += 1
                if first_mismatch_row is None:
                    first_mismatch_row = row_count
            gene, values = parse_first_values(line, spotcheck_cells)
            if gene in genes:
                raise RuntimeError(
                    "duplicate gene at row {0}".format(row_count)
                )
            genes.add(gene)
            for index, value in enumerate(values):
                observed_counts[index] += value
                if value > 0:
                    observed_features[index] += 1
            if row_count % 5000 == 0:
                print(
                    "GSE255460 stream audit: {0:,} gene rows".format(row_count)
                )
                sys.stdout.flush()

    result = {
        "member": args.member,
        "member_size": member.size,
        "n_cells": n_cells,
        "n_genes": row_count,
        "unique_genes": len(genes),
        "header_md5": header_md5,
        "integer_nonnegative_spotcheck": True,
        "observed_counts": observed_counts,
        "observed_features": observed_features,
        "structure_rows_checked": row_count,
        "structure_field_mismatches": field_mismatches,
        "first_structure_mismatch_row": first_mismatch_row,
        "expected_fields": expected_fields,
        "elapsed_seconds": time.time() - started,
    }
    temporary = args.output + ".tmp-{0}".format(os.getpid())
    with open(temporary, mode="w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if os.path.exists(args.output):
        os.remove(args.output)
    os.rename(temporary, args.output)
    print(
        "GSE255460 stream audit completed: {0:,} genes, {1:,} cells, "
        "{2:.1f} seconds".format(row_count, n_cells, result["elapsed_seconds"])
    )


if __name__ == "__main__":
    main()
