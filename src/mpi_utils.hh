#ifndef CLUSTRGAME_MPI_UTILS_H_
#define CLUSTRGAME_MPI_UTILS_H_


#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include "types.h"


static inline len_t get_numbefore(const len_local_t nrows_local, MPI_Comm comm)
{
  len_t nb4 = 0;
  int size, my_rank;
  
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &my_rank);
  
  for (int rank=1; rank<=size; rank++)
  {
    if (my_rank == (rank - 1))
    {
      len_t nb4_send = nb4 + nrows_local;
      MPI_Send(&nb4_send, 1, MPI_LEN_T, rank, 0, comm);
    }
    else if (my_rank == rank)
    {
      MPI_Status status;
      len_t nr_prev_rank;
      MPI_Recv(&nr_prev_rank, 1, MPI_LEN_T, rank-1, 0, comm, &status);
      
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
