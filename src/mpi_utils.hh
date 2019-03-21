#ifndef KMEANS_MPI_UTILS_H_
#define KMEANS_MPI_UTILS_H_


#ifdef __cplusplus
#define OMPI_SKIP_MPICXX 1
#endif
#include <mpi.h>


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



static inline int allreduce_real(const int len, float *const __restrict__ x, const MPI_Comm comm)
{
  return MPI_Allreduce(MPI_IN_PLACE, x, len, MPI_FLOAT, MPI_SUM, comm);
}

static inline int allreduce_real(const int len, double *const __restrict__ x, const MPI_Comm comm)
{
  return MPI_Allreduce(MPI_IN_PLACE, x, len, MPI_DOUBLE, MPI_SUM, comm);
}


#endif
