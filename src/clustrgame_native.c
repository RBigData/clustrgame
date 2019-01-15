/* Automatically generated. Do not edit by hand. */
  
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern SEXP R_clustrgame_init(SEXP comm_);
extern SEXP R_kmeans(SEXP data, SEXP m, SEXP n, SEXP k, SEXP maxiter, SEXP comm_);

static const R_CallMethodDef CallEntries[] = {
  {"R_clustrgame_init", (DL_FUNC) &R_clustrgame_init, 1},
  {"R_kmeans", (DL_FUNC) &R_kmeans, 6},
  {NULL, NULL, 0}
};

void R_init_clustrgame(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
