#include <mpi.h>
#include <Rinternals.h>

#include "../mpi_utils.h"


SEXP R_clustrgame_init(SEXP comm_)
{
  (void) comm_;
  return R_NilValue;
}
