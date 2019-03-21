#include <Rinternals.h>
#include <float/float32.h>

#include "cpu/kmeans_cpu.hh"

#include "common.h"
#include "kmeans_types.hh"
#include "Rmpi.h"


extern "C" SEXP R_kmeans(SEXP data, SEXP m, SEXP k_, SEXP maxiter, SEXP comm_)
{
  int check;
  SEXP ret, ret_names;
  SEXP ret_centers, ret_labels, ret_niters;
  
  MPI_Comm comm = get_mpi_comm_from_Robj(comm_);
  
  const int m_local = nrows(data);
  const int n = ncols(data);
  const int k = INTEGER(k_)[0];
  
  kmeans_opts opts;
  opts.k = k;
  opts.maxiter = INTEGER(maxiter)[0];
  opts.indexing = INDEXING_ONE;
  
  PROTECT(ret_labels = allocVector(INTSXP, m_local));
  PROTECT(ret_niters = allocVector(INTSXP, 1));
  
  if (TYPEOF(data) == REALSXP)
  {
    PROTECT(ret_centers = allocMatrix(REALSXP, n, k));
    
    shaq<double> x;
    x.nrows = INTEGER(m)[0];
    x.ncols = n;
    x.nrows_local = m_local;
    x.data = REAL(data);
    x.comm = comm;
    
    kmeans_vals<double> km;
    km.centers = REAL(ret_centers);
    km.labels = INTEGER(ret_labels);
    
    check = kmeans(&x, &km, &opts);
  }
  else if (TYPEOF(data) == INTSXP)
  {
    PROTECT(ret_centers = allocMatrix(INTSXP, n, k));
    
    shaq<float> x;
    x.nrows = INTEGER(m)[0];
    x.ncols = n;
    x.nrows_local = m_local;
    x.data = FLOAT(data);
    x.comm = comm;
    
    kmeans_vals<float> km;
    km.centers = FLOAT(ret_centers);
    km.labels = INTEGER(ret_labels);
    
    check = kmeans(&x, &km, &opts);
  }
  else
    error("this should be impossible; please contact the developers\n");
  
  if (check == ERROR_MPI)
    R_err_mpi(check, comm);
  else if (check == ERROR_MALLOC)
    R_err_malloc(comm);
  
  INTEGER(ret_niters)[0] = check;
  
  PROTECT(ret = allocVector(VECSXP, 3));
  PROTECT(ret_names = allocVector(STRSXP, 3));
  
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
