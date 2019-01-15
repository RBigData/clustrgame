#' @useDynLib clustrgame R_clustrgame_init
.onLoad = function(libname, pkgname)
{
  comm_ptr = pbdMPI::get.mpi.comm.ptr(.pbd_env$SPMD.CT$comm)
  .Call(R_clustrgame_init, comm_ptr)
  
  invisible()
}
