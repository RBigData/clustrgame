#ifndef KMEANS_CPU_H_
#define KMEANS_CPU_H_


#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <Rinternals.h>
// #include <float/float32.h>
// #include <float/slapack.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "../common.h"
#include "../restrict.h"
#include "../types.h"

#include "../kmeans.hh"
#include "../mpi_utils.hh"

#include "rand.h"
#include "utils.hh"


#define CHECKMPI(x) {if((x) != MPI_SUCCESS){return ERROR_MPI;}}
#define CHECKRET(x) {if((x) != ERROR_NONE){goto cleanup;}}


template <typename REAL>
static inline int kmeans_init(const shaq<REAL> *const restrict x, kmeans_vals<REAL> *const restrict km, const kmeans_opts *const restrict opts)
{
  const len_local_t m = NROWS_LOCAL(x);
  const len_local_t n = NCOLS(x);
  const int k = opts->k;
  
  int rank;
  MPI_Comm_rank(COMM(x), &rank);
  
  len_t *rows = (len_t*) malloc(k * sizeof(*rows));
  len_local_t *rows_local = (len_local_t*) malloc(k * sizeof(*rows_local));
  if (rows == NULL || rows_local == NULL)
  {
    FREE(rows);
    FREE(rows_local);
    return ERROR_MALLOC;
  }
  
  uint64_t nb4 = get_numbefore(m, x->comm);
  
  reservoir_sampler(NROWS(x), k, rows);
  sort_insertion(k, rows);
  
  for (int i=0; i<k; i++)
    rows_local[i] = nb4 <= rows[i] && rows[i] < nb4+m;
  
  for (int row=0; row<k; row++)
  {
    if (rows_local[row])
    {
      for (len_local_t j=0; j<n; j++)
        km->centers[j + row*n] = DATA(x)[(rows[row]-nb4) + m*j];
    }
  }
  
  free(rows);
  free(rows_local);
  
  int check = allreduce_real(n*k, km->centers, COMM(x));
  CHECKMPI(check);
  return ERROR_NONE;
}



template <typename REAL>
static inline int kmeans_update(const shaq<REAL> *const restrict x, kmeans_vals<REAL> *const restrict km, const kmeans_opts *const restrict opts)
{
  int check;
  const len_local_t m = NROWS_LOCAL(x);
  const len_local_t n = NCOLS(x);
  const int k = opts->k;
  
  REAL *const restrict centers = km->centers;
  int *const restrict labels = km->labels;
  int *const restrict nlabels = km->nlabels;
  
  memset(centers, 0, n*k*sizeof(*centers));
  for (len_local_t j=0; j<n; j++)
  {
    for (len_local_t i=0; i<m; i++)
      centers[j + n*labels[i]] += DATA(x)[i + m*j];
  }
  
  check = allreduce_real(n*k, centers, COMM(x));
  CHECKMPI(check);
  
  
  memset(nlabels, 0, k*sizeof(*nlabels));
  for (int i=0; i<m; i++)
    nlabels[labels[i]]++;
  
  check = MPI_Allreduce(MPI_IN_PLACE, nlabels, k, MPI_INTEGER, MPI_SUM, COMM(x));
  CHECKMPI(check);
  
  for (int j=0; j<k; j++)
  {
    for (len_local_t i=0; i<n; i++)
      centers[i + n*j] /= (double)nlabels[j];
  }
  
  return ERROR_NONE;
}



template <typename REAL>
static inline int kmeans_assign_single(const int m, const int n, const int k, const REAL *const restrict x, const REAL *const restrict centers)
{
  REAL min = INFINITY;
  int min_ind = -1;
  
  for (int j=0; j<k; j++)
  {
    REAL test = 0.0;
    
    for (len_local_t i=0; i<n; i++)
    {
      const REAL tmp = x[m*i] - centers[i + n*j];
      test += tmp*tmp;
    }
    
    if (test < min || min_ind == -1)
    {
      min = test;
      min_ind = j;
    }
  }
  
  return min_ind;
}



template <typename REAL>
static inline void kmeans_assign(const shaq<REAL> *const restrict x, kmeans_vals<REAL> *const restrict km, const kmeans_opts *const restrict opts)
{
  const len_local_t m = NROWS_LOCAL(x);
  const len_local_t n = NCOLS(x);
  const int k = opts->k;
  
  for (len_local_t i=0; i<m; i++)
    km->labels[i] = kmeans_assign_single(m, n, k, DATA(x)+i, km->centers);
}



// returns number of iterations
template <typename REAL>
static inline int kmeans(const shaq<REAL> *const restrict x, kmeans_vals<REAL> *const restrict km, const kmeans_opts *const restrict opts)
{
  int ret;
  int niters;
  const int k = opts->k;
  const len_local_t nk = NCOLS(x) * k;
  
  int *nlabels = (int*) malloc(k * sizeof(*nlabels));
  km->nlabels = nlabels;
  
  REAL *centers_old = (REAL*) malloc(nk * sizeof(*centers_old));
  memset(centers_old, 0, nk*sizeof(*centers_old));
  
  if (nlabels == NULL || centers_old == NULL)
  {
    ret = ERROR_MALLOC;
    goto cleanup;
  }
  
  ret = kmeans_init(x, km, opts);
  CHECKRET(ret);
  kmeans_assign(x, km, opts);
  
  for (niters=0; niters<opts->maxiter; niters++)
  {
    ret = kmeans_update(x, km, opts);
    CHECKRET(ret);
    kmeans_assign(x, km, opts);
    
    int check = all_equal(nk, km->centers, centers_old);
    if (check == TRUE)
      break;
    
    memcpy(centers_old, km->centers, nk*sizeof(*centers_old));
  }
  
  
  if (!opts->zero_indexed)
    add1(NROWS_LOCAL(x), km->labels);
  
  if (niters < opts->maxiter)
    niters++;
  
  ret = niters;
  
  
cleanup:
  FREE(centers_old);
  FREE(nlabels);
  km->nlabels = NULL;
  
  return ret;
}


#endif
