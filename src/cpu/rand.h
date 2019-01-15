#ifndef KMEANS_RAND_H_
#define KMEANS_RAND_H_


#include <R_ext/Random.h>
#include "types.h"


static inline int unif_rand_int(const int low, const int high)
{
  return (int) low + (high + 1 - low)*unif_rand();
}



static inline int res_sampler(const len_t nlines_in, const len_t nlines_out, len_t *samp)
{
  len_t i, j;
  
  for (i=0; i<nlines_out; i++)
    samp[i] = i + 1;
  
  GetRNGstate();
  
  for (i=nlines_out; i<nlines_in; i++)
  {
    j = unif_rand_int(0, i-1);
    if (j < nlines_out)
      samp[j] = i + 1;
  }
  
  PutRNGstate();
  
  return 0;
}


#endif
