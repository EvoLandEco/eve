job_status_parser <- function(file = NULL) {
  raw_log <- readChar(file, file.info(file)$size)

  set_id <-
    gsub(".+(Starting parameter set )(\\d+)(.+)", "\\2", raw_log)
  status <-
    gsub(".+(State\\s+\\:\\s)(.+)\n(Submit).+", "\\2", raw_log)

  if (nchar(status) >= 30 | nchar(status) < 1) {
    set_id <- "Unknown"
    status <- "RUNNING"
  }

  return(c(set_id, status))
}


#' @name job_status
#' @title Parsing log files in a folder into a matrix of job status
#' @description Function to parse outputs in all log files of a folder, then
#' return a matrix to show the status of the jobs
#' @param path path to the result folder
#' @param which specify status to be filtered
#' @author Tianjian Qin
#' @export job_status
job_status <- function(path = NULL, which = NULL) {
  if (dir.exists(path)) {
    file_names <- list.files(file.path(path, "logs"))
    status <-
      lapply(X = file.path(path, "logs", file_names),
             FUN = job_status_parser)
  } else {
    stop("Invalid path")
  }

  status <- matrix(unlist(status), ncol = 2, byrow = TRUE)
  status[, 1] <- seq.int(nrow(status))
  summary <- table(status[, 2])

  if (!is.null(which)) {
    status <- status[status[, 2] == which, ]
  }

  if (is.null(nrow(status))) {
    names(status) <- c("Set", "Status")
  } else if (nrow(status) == 0) {
    message("No job with status ", which)
  } else {
    colnames(status) <- c("Set", "Status")
  }

  return(list(status = status,
              summary = summary))
}



job_detail_parser <- function(file = NULL) {
  raw_string <- readChar(file, file.info(file)$size)
  match_pattern <-
    "(Job ID)\\s+:\\s(\\d+)\n(Name)\\s+:\\s(.+)\n(User)\\s+:\\s(.+)\n(Partition)\\s+:\\s(.+)\n(Nodes)\\s+:\\s(.+)\n(Number of Nodes)\\s+:\\s(.+)\n(Cores)\\s+:\\s(.+)\n(Number of Tasks)\\s+:\\s(.+)\n(State)\\s+:\\s(.+)\n(Submit)\\s+:\\s(.+)\n(Start)\\s+:\\s(.+)\n(End)\\s+:\\s(.+)\n(Reserved walltime)\\s+:\\s(.+)\n(Used walltime)\\s+:\\s(.+)\n(Used CPU time)\\s+:\\s(.+)\n(% User \\(Computation\\))\\s*:\\s(.+)\n(% System \\(I\\/O\\))\\s+:\\s(.+)\n(Mem reserved)\\s+:\\s(.+)\n(Max Mem \\(Node\\/step\\))\\s+:\\s(.+)\n(Full Max Mem usage)\\s+:\\s(.+)\n(Total Disk Read)\\s+:\\s(.+)\n(Total Disk Write)\\s+:\\s(.+)\n"
  matched_strings <- stringr::str_match(string = raw_string,
                                        pattern = match_pattern)
  detail <-
    matrix(data = matched_strings[-1],
           ncol = 2,
           byrow = TRUE)

  if (NA %in% detail) {
    warning(paste0("Failed to parse ", file, "\nIs the job already finished?"))
    return(NA)
  } else {
    return(detail)
  }
}


#' @name job_detail
#' @title Parsing log files in a folder into a list of job details
#' @description Function to parse job details in all log files of a folder, then
#' return a list storing those information
#' @param path path to the result folder
#' @param which specify which job to be parsed
#' @author Tianjian Qin
#' @export job_detail
job_detail <- function(path = NULL, which = NULL) {
  if (dir.exists(path)) {
    file_names <- list.files(file.path(path, "logs"))
    detail <-
      lapply(X = file.path(path, "logs",file_names),
             FUN = job_detail_parser)
  } else {
    stop("Invalid path")
  }

  if (is.null(which)) {
      return(detail)
  } else {
      return(detail[[which]])
  }
}



#' @name job_parameter
#' @title Find out which parameter set is used for a job
#' @description Function to find out which parameter set is used for a job
#' @param path path to the result folder
#' @param which specify which job to be parsed
#' @author Tianjian Qin
#' @export job_detail
job_parameter <- function(path = NULL, which = NULL) {
  params <- read.table(file.path(path, "params"))

  if (!is.null(which)) {
    return(params[which, ])
  } else {
    return(params)
  }
}