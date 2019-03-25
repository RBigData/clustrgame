#' @useDynLib clustrgame R_clustrgame_init
.onLoad = function(libname, pkgname)
{
  s = search()
  if ("package:dimrgame" %in% s || "package:glmrgame" %in% s)
    return(invisible())
  
  comm_ptr = pbdMPI::get.mpi.comm.ptr(.pbd_env$SPMD.CT$comm)
  .Call(R_clustrgame_init, comm_ptr)
  
  invisible()
}
