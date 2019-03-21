#ifndef CLUSTRGAME_KMEANS_TYPES_H_
#define CLUSTRGAME_KMEANS_TYPES_H_


#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <cstdint>

#include "types.h"

#define INDEXING_ZERO 0
#define INDEXING_ONE 1

#define NROWS(x) (x->nrows)
#define NCOLS(x) (x->ncols)
#define NROWS_LOCAL(x) (x->nrows_local)
#define NCOLS_LOCAL(x) (x->ncols)
#define DATA(x) (x->data)
#define COMM(x) (x->comm)

template <typename REAL>
struct shaq
{
  len_t nrows;
  len_local_t ncols;
  len_local_t nrows_local;
  REAL *restrict data;
  MPI_Comm comm;
};

template <typename REAL>
struct kmeans_vals
{
  REAL *restrict centers;
  int *restrict labels;
  int *restrict nlabels;
};

typedef struct kmeans_opts
{
  int k;
  int maxiter;
  int seed;
  int indexing;
} kmeans_opts;


#endif
