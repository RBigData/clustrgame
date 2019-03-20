#include <Rinternals.h>
// #include <float/float32.h>
// #include <float/slapack.h>

#include "cpu/kmeans_cpu.hh"

#include "common.h"
#include "kmeans.hh"


extern "C" SEXP R_kmeans(SEXP data, SEXP m, SEXP k_, SEXP maxiter, SEXP comm_)
{
  SEXP ret, ret_names;
  SEXP ret_centers, ret_labels, ret_niters;
  shaq<double> x;
  kmeans_vals<double> km;
  kmeans_opts opts;
  
  MPI_Comm comm = get_mpi_comm_from_Robj(comm_);
  
  const int m_local = nrows(data);
  const int n = ncols(data);
  const int k = INTEGER(k_)[0];
  
  PROTECT(ret_centers = allocMatrix(REALSXP, n, k));
  PROTECT(ret_labels = allocVector(INTSXP, m_local));
  PROTECT(ret_niters = allocVector(INTSXP, 1));
  
  PROTECT(ret = allocVector(VECSXP, 3));
  PROTECT(ret_names = allocVector(STRSXP, 3));
  
  x.nrows = INTEGER(m)[0];
  x.ncols = n;
  x.nrows_local = m_local;
  x.data = REAL(data);
  x.comm = comm;
  
  km.centers = REAL(ret_centers);
  km.labels = INTEGER(ret_labels);
  
  opts.k = k;
  opts.maxiter = INTEGER(maxiter)[0];
  opts.zero_indexed = 0;
  
  int check = kmeans(&x, &km, &opts);
  
  if (check == ERROR_MPI)
    throw_err_mpi(check, comm);
  else if (check == ERROR_MALLOC)
    throw_err_malloc(comm);
  
  INTEGER(ret_niters)[0] = check;
  
  SET_VECTOR_ELT(ret, 0, ret_centers);
  SET_VECTOR_ELT(ret, 1, ret_labels);
  SET_VECTOR_ELT(ret, 2, ret_niters);
  SET_STRING_ELT(ret_names, 0, mkChar("centers"));
  SET_STRING_ELT(ret_names, 1, mkChar("labels"));
  SET_STRING_ELT(ret_names, 2, mkChar("niters"));
  setAttrib(ret, R_NamesSymbol, ret_names);
  
  UNPROTECT(5);
  return ret;
}
