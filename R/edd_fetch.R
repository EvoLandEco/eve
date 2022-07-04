#' Function to fetch simulation files from an \code{FTP} or \code{SFTP} server.
#'
#' @param url \code{URL}s of the files to be downloaded.
#' @param credentials Credentials for a \code{FTP} or \code{SFTP} server. Do not
#' use \code{credentials} if the server does not require authentication.
#' \code{credentials} takes the format: \code{"username:password"}.
#' @author Tianjian Qin
#' @export
edd_fetch <- function(url = NULL, credentials = NULL) {
  edd_fetch_logs(url = NULL, credentials = NULL)
  edd_fecth_results(url = NULL, credentials = NULL)
}



#' Function to fetch log files from an \code{FTP} or \code{SFTP} server.
#'
#' @param url \code{URL}s of the files to be downloaded.
#' @param credentials Credentials for a \code{FTP} or \code{SFTP} server. Do not
#' use \code{credentials} if the server does not require authentication.
#' \code{credentials} takes the format: \code{"username:password"}.
#' @author Tianjian Qin
#' @export
edd_fetch_logs <- function(url = NULL, credentials = NULL) {

}



#' Function to fetch simulation results from an \code{FTP} or \code{SFTP} server.
#'
#' @param url \code{URL}s of the files to be downloaded.
#' @param credentials Credentials for a \code{FTP} or \code{SFTP} server. Do not
#' use \code{credentials} if the server does not require authentication.
#' \code{credentials} takes the format: \code{"username:password"}.
#' @author Tianjian Qin
#' @export
edd_fetch_results <- function(url = NULL, credentials = NULL) {

}



