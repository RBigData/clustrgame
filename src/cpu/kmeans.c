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


#define MALLOC_FAIL -1



static inline void add1(const int n, int *const x)
{
  for (int i=0; i<n; i++)
    x[i]++;
}



static inline void kmeans_init(const shaq *const restrict x, kmeans_vals *const restrict km, const kmeans_opts *const restrict opts)
{
  const int n = NROWS_LOCAL(x);
  const int k = opts->k;
  
  int rank;
  MPI_Comm_rank(COMM(x), &rank);
  
  len_t *rows = malloc(k * sizeof(*rows));
  int *rows_local = malloc(k * sizeof(*rows_local));
  
  uint64_t nb4 = get_numbefore(NROWS_LOCAL(x), x->comm);
  
  reservoir_sampler(NROWS(x), k, rows);
  sort_insertion(k, rows);
  
  for (int i=0; i<k; i++)
    rows_local[i] = nb4 <= rows[i] && rows[i] < nb4+NROWS_LOCAL(x);
  
  for (int row=0; row<k; row++)
  {
    if (rows_local[row])
    {
      for (int j=0; j<NCOLS(x); j++)
        km->centers[j + row*NCOLS(x)] = DATA(x)[(rows[row]-nb4) + NROWS_LOCAL(x)*j];
    }
  }
  
  free(rows);
  free(rows_local);
  
  MPI_Barrier(COMM(x));
  
  int check = MPI_Allreduce(MPI_IN_PLACE, km->centers, n*k, MPI_DOUBLE, MPI_SUM, COMM(x));
  // TODO
  // if (check != MPI_SUCCESS)
}



static inline void kmeans_update(const shaq *const restrict x, kmeans_vals *const restrict km, const kmeans_opts *const restrict opts)
{
  int check;
  const int m = NROWS_LOCAL(x);
  const int n = NCOLS_LOCAL(x);
  const int k = opts->k;
  
  memset(km->centers, 0, n*k*sizeof(*km->centers));
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<m; i++)
      km->centers[j + n*km->labels[i]] += DATA(x)[i + m*j];
  }
  
  check = MPI_Allreduce(MPI_IN_PLACE, km->centers, n*k, MPI_DOUBLE, MPI_SUM, x->comm);
  if (check != MPI_SUCCESS)
  {
    free(km->nlabels);
    R_mpi_throw_err(check, x->comm);
  }
  
  
  memset(km->nlabels, 0, k*sizeof(*km->nlabels));
  for (int i=0; i<m; i++)
    km->nlabels[km->labels[i]]++;
  
  check = MPI_Allreduce(MPI_IN_PLACE, km->nlabels, k, MPI_INTEGER, MPI_SUM, x->comm);
  if (check != MPI_SUCCESS)
  {
    free(km->nlabels);
    R_mpi_throw_err(check, x->comm);
  }
    
  for (int j=0; j<k; j++)
  {
    for (int i=0; i<n; i++)
      km->centers[i + n*j] /= (double)km->nlabels[j];
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




static inline int kmeans_assign(const shaq *const restrict x, kmeans_vals *const restrict km, const kmeans_opts *const restrict opts)
{
  const int m = NROWS_LOCAL(x);
  const int n = NCOLS_LOCAL(x);
  const int k = opts->k;
  
  for (int i=0; i<m; i++)
    km->labels[i] = kmeans_assign_single(m, n, k, x->data+i, km->centers);
}



// returns number of iterations
static inline int kmeans(const shaq *const restrict x, kmeans_vals *const restrict km, const kmeans_opts *const restrict opts)
{
  int niters;
  kmeans_init(x, km, opts);
  
  exit(1);
  
  int *nlabels = malloc(opts->k * sizeof(*nlabels));
  if (nlabels == NULL)
    error("OOM");
  km->nlabels = nlabels;
  
  
  for (niters=0; niters<opts->maxiter; niters++)
  {
    kmeans_update(x, km, opts);
    
    // TODO check for convergence
  }
  
  
  free(nlabels);
  km->nlabels = NULL;
  
  if (!opts->zero_indexed)
    add1(x->nrows_local, km->labels);
  
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
