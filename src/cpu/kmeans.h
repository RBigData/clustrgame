#ifndef KMEANS_H_
#define KMEANS_H_


#include <mpi.h>
#include <stdint.h>

#include "types.h"

#define NROWS(x) (x->nrows)
#define NCOLS(x) (x->ncols)
#define NROWS_LOCAL(x) (x->nrows_local)
#define NCOLS_LOCAL(x) (x->ncols)
#define DATA(x) (x->data)
#define COMM(x) (x->comm)

typedef struct shaq
{
  len_t nrows;
  int ncols;
  int nrows_local;
  double *restrict data;
  MPI_Comm comm;
} shaq;

typedef struct kmeans_vals
{
  double *restrict centers;
  int *restrict labels;
  int *restrict nlabels;
} kmeans_vals;

typedef struct kmeans_opts
{
  int k;
  int maxiter;
  int seed;
  int zero_indexed;
} kmeans_opts;


#endif
