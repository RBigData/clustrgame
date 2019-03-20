#include <cuda_runtime.h>
#include <mpi.h>

extern "C" {
#include <Rinternals.h>
#include "../mpi_utils.h"
}


extern "C" SEXP R_clustrgame_init(SEXP comm_)
{
  int ngpus;
  int rank;
  int id;
  
  MPI_Comm comm = get_mpi_comm_from_Robj(comm_);
  MPI_Comm_rank(comm, &rank);
  
  cudaGetDeviceCount(&ngpus);
  
  id = rank % ngpus;
  cudaSetDevice(id);
  
  return R_NilValue;
}
