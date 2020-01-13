.onUnload <- function (libpath) {
  library.dynam.unload("seagull", libpath)
}