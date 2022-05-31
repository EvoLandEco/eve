#' save_result
#'
#' @author Tianjian Qin
save_result <- function(out,
                        name,
                        set = NULL) {
  if (is.null(set)) {
    file_name <- paste0(name, ".RData")
  } else {
    file_name <- paste0(name, "_", set, ".RData")
  }

  out_folder <- file.path(getwd(),
                          "result",
                          name)

  file_path <- file.path(out_folder, file_name)

  message(paste0("\nSaving result to ",
                 file_path,
                 "\n"))

  save(out, file = file_path)

  if (file.exists(file_path)) {
    message(paste0("Result saved", "\n"))
  } else {
    warning(paste0("Failed to save result to ",
                   file_path))
  }
}



#' check_folder
#'
#' @author Tianjian Qin
check_folder <- function(name, verbose = TRUE) {
  out_folder <- file.path(getwd(),
                          "result",
                          name)
  if (!dir.exists(out_folder)) {
    if (verbose == TRUE) {
      message(paste0("Attempting to create folder ",
                     out_folder,
                     "\n"))
    }
    dir.create(out_folder, recursive = TRUE)
  } else {
    if (verbose == TRUE) {
      message(out_folder, " already exists\n")
    }
  }

  if (!dir.exists(out_folder)) {
    stop(paste0("Failed to create folder", out_folder, "\n"))
  } else {
    if (verbose == TRUE) {
      message("Folder created\n")
    }
  }
}
