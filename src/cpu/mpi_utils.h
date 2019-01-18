#ifndef KMEANS_MPI_UTILS_H_
#define KMEANS_MPI_UTILS_H_


#include <mpi.h>


#define MPI_CHECK(comm, check) if (check != MPI_SUCCESS) R_mpi_throw_err(check, comm);

static inline void throw_err_mpi(int check, const MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0)
    error("MPI_Allreduce returned error code %d\n", check);
  else
    error(""); // FIXME
}



static inline void throw_err_malloc(const MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0)
    error("Out of memory\n");
  else
    error(""); // FIXME
}



static inline MPI_Comm* get_mpi_comm_from_Robj(SEXP comm_)
{
  MPI_Comm *comm = (MPI_Comm*) R_ExternalPtrAddr(comm_);
  return comm;
}



static inline uint64_t get_numbefore(const int nrows_local, MPI_Comm comm)
{
  uint64_t nb4 = 0;
  int size, my_rank;
  
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &my_rank);
  
  for (int rank=1; rank<=size; rank++)
  {
    if (my_rank == (rank - 1))
    {
      uint64_t nb4_send = nb4 + nrows_local;
      MPI_Send(&nb4_send, 1, MPI_UINT64_T, rank, 0, comm);
    }
    else if (my_rank == rank)
    {
      MPI_Status status;
      uint64_t nr_prev_rank;
      MPI_Recv(&nr_prev_rank, 1, MPI_UINT64_T, rank-1, 0, comm, &status);
      
      nb4 += nr_prev_rank;
    }
  }
  
  return nb4;
}


#endif
