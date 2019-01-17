#' km
#' 
#' k-means via Lloyd's Algorithm.
#' 
#' Note that the function does not respect \code{set.seed()} or
#' \code{comm.set.seed()}.  For managing random seeds, use the \code{seed}
#' parameter.
#' 
#' @details
#' The iterations stop either when the maximum number of iterations have been
#' achieved, or when the centers in the current iteration are basically the same
#' (within \code{1e-8}) as the centers from the previous iteration.
#' 
#' For best performance, the data should be as balanced as possible across all
#' MPI ranks.
#' 
#' @section Communication:
#' Most of the computation is local. However, at each iteration there is a
#' length \code{n*k} and a length \code{k} allreduce call to update the centers.
#' There is also a check at the beginning of the call to find out how many
#' observations come before the current process's data, which is an allgather
#' operation.
#' 
#' @param x
#' A shaq.
#' @param k
#' The 'k' in k-means.
#' @param maxiter
#' The maximum number of iterations possible.
#' @param seed
#' A seed for determining the (random) initial centroids.  Each process has to
#' use the same seed or very strange things may happen.  If you do not provide
#' a seed, a good initial seed will be chosen.
#' 
#' @return
#' A list containing the cluster centers (global), the observation labels i.e.
#' the assignments to clusters (distributed shaq), and the total number of
#' iterations (global).
#'
#' @references
#' Phillips, J.. Data Mining: Algorithms, Geometry, and Probability.
#' https://www.cs.utah.edu/~jeffp/DMBook/DM-AGP.html
#' 
#' @useDynLib clustrgame R_kmeans
#' @export
km_game = function(x, k=2, maxiter=100, seed=1234)
{
  kazaam:::check.is.shaq(x)
  kazaam:::check.is.posint(k)
  kazaam:::check.is.posint(maxiter)
  kazaam:::check.is.posint(seed)
  
  k = as.integer(k)
  maxiter = as.integer(maxiter)
  comm_ptr = pbdMPI::get.mpi.comm.ptr(.pbd_env$SPMD.CT$comm)
  
  data = DATA(x)
  if (!is.double(data))
    storage.mode(data) = "double"
  
  m = nrow(x)
  m = as.integer(m) # FIXME
  
  ret = .Call(R_kmeans, data, m, k, maxiter, comm_ptr)
  ret
}
