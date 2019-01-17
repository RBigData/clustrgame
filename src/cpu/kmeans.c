// #include <float/float32.h>
// #include <float/slapack.h>
#include <mpi.h>
#include <Rinternals.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "kmeans.h"
#include "mpi_utils.h"
#include "rand.h"
#include "utils.h"


#define MALLOC_FAIL -1


static inline void kmeans_init(const shaq *const restrict x, kmeans_vals *const restrict km, const kmeans_opts *const restrict opts)
{
  const int m = NROWS_LOCAL(x);
  const int n = NCOLS(x);
  const int k = opts->k;
  
  int rank;
  MPI_Comm_rank(COMM(x), &rank);
  
  len_t *rows = malloc(k * sizeof(*rows));
  int *rows_local = malloc(k * sizeof(*rows_local));
  
  uint64_t nb4 = get_numbefore(m, x->comm);
  
  reservoir_sampler(NROWS(x), k, rows);
  sort_insertion(k, rows);
  
  for (int i=0; i<k; i++)
    rows_local[i] = nb4 <= rows[i] && rows[i] < nb4+m;
  
  for (int row=0; row<k; row++)
  {
    if (rows_local[row])
    {
      for (int j=0; j<n; j++)
        km->centers[j + row*n] = DATA(x)[(rows[row]-nb4) + m*j];
    }
  }
  
  free(rows);
  free(rows_local);
  
  int check = MPI_Allreduce(MPI_IN_PLACE, km->centers, n*k, MPI_DOUBLE, MPI_SUM, COMM(x));
  // TODO
  // if (check != MPI_SUCCESS)
}



static inline void kmeans_update(const shaq *const restrict x, kmeans_vals *const restrict km, const kmeans_opts *const restrict opts)
{
  int check;
  const int m = NROWS_LOCAL(x);
  const int n = NCOLS(x);
  const int k = opts->k;
  
  double *const restrict centers = km->centers;
  int *const restrict labels = km->labels;
  int *const restrict nlabels = km->nlabels;
  
  memset(centers, 0, n*k*sizeof(*centers));
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<m; i++)
      centers[j + n*labels[i]] += DATA(x)[i + m*j];
  }
  
  check = MPI_Allreduce(MPI_IN_PLACE, centers, n*k, MPI_DOUBLE, MPI_SUM, COMM(x));
  // if (check != MPI_SUCCESS)
  // {
  //   free(nlabels);
  //   R_mpi_throw_err(check, COMM(x));
  // }
  
  
  memset(nlabels, 0, k*sizeof(*nlabels));
  for (int i=0; i<m; i++)
    nlabels[labels[i]]++;
  
  check = MPI_Allreduce(MPI_IN_PLACE, nlabels, k, MPI_INTEGER, MPI_SUM, COMM(x));
  // if (check != MPI_SUCCESS)
  // {
  //   free(nlabels);
  //   R_mpi_throw_err(check, COMM(x));
  // }
    
  for (int j=0; j<k; j++)
  {
    for (int i=0; i<n; i++)
      centers[i + n*j] /= (double)nlabels[j];
  }
}



static inline int kmeans_assign_single(const int m, const int n, const int k, const double *const restrict x, const double *const restrict centers)
{
  double min = INFINITY;
  int min_ind = -1;
  
  for (int j=0; j<k; j++)
  {
    double test = 0.0;
    
    for (int i=0; i<n; i++)
    {
      const double tmp = x[m*i] - centers[i + n*j];
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



static inline void kmeans_assign(const shaq *const restrict x, kmeans_vals *const restrict km, const kmeans_opts *const restrict opts)
{
  const int m = NROWS_LOCAL(x);
  const int n = NCOLS(x);
  const int k = opts->k;
  
  for (int i=0; i<m; i++)
    km->labels[i] = kmeans_assign_single(m, n, k, DATA(x)+i, km->centers);
}



// returns number of iterations
static inline int kmeans(const shaq *const restrict x, kmeans_vals *const restrict km, const kmeans_opts *const restrict opts)
{
  int niters;
  const int k = opts->k;
  const int nk = NCOLS(x) * k;
  
  int *nlabels = malloc(k * sizeof(*nlabels));
  if (nlabels == NULL)
    error("OOM");
  km->nlabels = nlabels;
  
  double *centers_old = malloc(nk * sizeof(*centers_old));
  memset(centers_old, 0, nk*sizeof(*centers_old));
  
  kmeans_init(x, km, opts);
  kmeans_assign(x, km, opts);
  
  for (niters=0; niters<opts->maxiter; niters++)
  {
    kmeans_update(x, km, opts);
    kmeans_assign(x, km, opts);
    
    int check = all_equal(nk, km->centers, centers_old);
    if (check == TRUE)
      break;
    
    memcpy(centers_old, km->centers, nk*sizeof(*centers_old));
  }
  
  
  free(centers_old);
  free(nlabels);
  km->nlabels = NULL;
  
  if (!opts->zero_indexed)
    add1(NROWS_LOCAL(x), km->labels);
  
  if (niters < opts->maxiter)
    niters++;
  
  return niters;
}



SEXP R_kmeans(SEXP data, SEXP m, SEXP k_, SEXP maxiter, SEXP comm_)
{
  SEXP ret, ret_names;
  SEXP centers, labels, niters;
  shaq x;
  kmeans_vals km;
  kmeans_opts opts;
  
  MPI_Comm comm = *(get_mpi_comm_from_Robj(comm_));
  
  const int m_local = nrows(data);
  const int n = ncols(data);
  const int k = INTEGER(k_)[0];
  
  PROTECT(centers = allocMatrix(REALSXP, n, k));
  PROTECT(labels = allocVector(INTSXP, m_local));
  PROTECT(niters = allocVector(INTSXP, 1));
  
  PROTECT(ret = allocVector(VECSXP, 3));
  PROTECT(ret_names = allocVector(STRSXP, 3));
  
  x.nrows = INTEGER(m)[0];
  x.ncols = n;
  x.nrows_local = m_local;
  x.data = REAL(data);
  x.comm = comm;
  
  km.centers = REAL(centers);
  km.labels = INTEGER(labels);
  
  opts.k = k;
  opts.maxiter = INTEGER(maxiter)[0];
  opts.zero_indexed = 0;
  
  INTEGER(niters)[0] = kmeans(&x, &km, &opts);
  
  SET_VECTOR_ELT(ret, 0, centers);
  SET_VECTOR_ELT(ret, 1, labels);
  SET_VECTOR_ELT(ret, 2, niters);
  SET_STRING_ELT(ret_names, 0, mkChar("centers"));
  SET_STRING_ELT(ret_names, 1, mkChar("labels"));
  SET_STRING_ELT(ret_names, 2, mkChar("niters"));
  setAttrib(ret, R_NamesSymbol, ret_names);
  
  UNPROTECT(5);
  return ret;
}
