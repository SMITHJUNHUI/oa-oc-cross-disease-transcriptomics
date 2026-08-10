`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

parse_cli_args <- function(args) {
  result <- list(
    mode = "full",
    config = "config/config.yml",
    force = FALSE,
    from = "",
    to = ""
  )

  for (arg in args) {
    if (identical(arg, "--force")) {
      result$force <- TRUE
    } else if (grepl("^--mode=", arg)) {
      result$mode <- sub("^--mode=", "", arg)
    } else if (grepl("^--config=", arg)) {
      result$config <- sub("^--config=", "", arg)
    } else if (grepl("^--from=", arg)) {
      result$from <- sub("^--from=", "", arg)
    } else if (grepl("^--to=", arg)) {
      result$to <- sub("^--to=", "", arg)
    } else {
      stop("Unknown command-line argument: ", arg, call. = FALSE)
    }
  }

  allowed_modes <- c("preflight", "core", "full", "tests")
  if (!result$mode %in% allowed_modes) {
    stop(
      "Invalid mode '", result$mode,
      "'. Choose one of: ", paste(allowed_modes, collapse = ", "),
      call. = FALSE
    )
  }
  result
}

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]", path) || startsWith(path, "/") ||
    startsWith(path, "\\\\")
}

resolve_path <- function(path, project_root, must_work = FALSE) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) {
    return(path)
  }

  expanded <- path.expand(path)
  candidate <- if (is_absolute_path(expanded)) {
    expanded
  } else {
    file.path(project_root, expanded)
  }

  normalizePath(
    candidate,
    winslash = "/",
    mustWork = must_work
  )
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(path)) {
      stop("Could not create directory: ", path, call. = FALSE)
    }
  }
  invisible(path)
}

initialize_logging <- function(output_dir) {
  log_dir <- ensure_dir(file.path(output_dir, "logs"))
  log_file <- file.path(
    log_dir,
    paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
  )
  options(ocoa.log_file = log_file)
  writeLines(character(), log_file, useBytes = TRUE)
  log_file
}

log_message <- function(level, ..., .sep = "") {
  text <- paste(..., sep = .sep, collapse = "")
  line <- sprintf(
    "[%s] %-5s %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    toupper(level),
    text
  )
  cat(line, "\n")

  log_file <- getOption("ocoa.log_file", NULL)
  if (!is.null(log_file)) {
    connection <- file(log_file, open = "a", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    writeLines(line, connection, useBytes = TRUE)
  }
  invisible(line)
}

log_info <- function(...) log_message("INFO", ...)
log_warn <- function(...) log_message("WARN", ...)
log_error <- function(...) log_message("ERROR", ...)

require_namespace <- function(package, reason = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    suffix <- if (is.null(reason)) "" else paste0(" (", reason, ")")
    stop(
      "Required R package is not installed: ", package, suffix,
      ". Run: Rscript setup.R --install",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

atomic_save_rds <- function(object, path, compress = TRUE) {
  ensure_dir(dirname(path))
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (file.exists(path)) {
    unlink(path, force = TRUE)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not atomically write checkpoint: ", path, call. = FALSE)
  }
  invisible(path)
}

write_utf8 <- function(lines, path) {
  ensure_dir(dirname(path))
  connection <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  writeLines(enc2utf8(lines), connection, useBytes = TRUE)
  invisible(path)
}

clean_filename <- function(value) {
  cleaned <- gsub("[^A-Za-z0-9._-]+", "_", value)
  gsub("^_+|_+$", "", cleaned)
}

human_bytes <- function(bytes) {
  units <- c("B", "KB", "MB", "GB", "TB")
  if (is.na(bytes) || bytes < 0) return(NA_character_)
  index <- min(floor(log(max(bytes, 1), 1024)) + 1L, length(units))
  sprintf("%.2f %s", bytes / 1024^(index - 1L), units[[index]])
}

safe_write_csv <- function(x, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  utils::write.csv(
    x,
    file = path,
    row.names = row.names,
    fileEncoding = "UTF-8",
    na = ""
  )
  invisible(path)
}

safe_numeric_matrix <- function(x, context = "matrix") {
  matrix <- as.matrix(x)
  storage.mode(matrix) <- "double"
  if (any(is.infinite(matrix), na.rm = TRUE)) {
    stop(context, " contains infinite values.", call. = FALSE)
  }
  matrix
}

compact_error <- function(error) {
  paste(class(error)[[1L]], conditionMessage(error), sep = ": ")
}
