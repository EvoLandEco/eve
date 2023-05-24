# Make cached versions of those computational-expensive functions upon package loading
.onLoad <- function(libname, pkgname) {
  edd_stat_cached <<- memoise::memoise(edd_stat)
}