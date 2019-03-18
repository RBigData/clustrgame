#ifndef KMEANS_RAND_H_
#define KMEANS_RAND_H_


#include <R_ext/Random.h>

#include "../types.h"


static inline int unif_rand_int(const int low, const int high)
{
  return (int) low + (high + 1 - low)*unif_rand();
}



#define INDEX_BASE 0
static inline int reservoir_sampler(const len_t m, const len_t k, len_t *samp)
{
  len_t i, j;
  
  for (i=0; i<k; i++)
    samp[i] = i + INDEX_BASE;
  
  GetRNGstate();
  
  for (i=k; i<m; i++)
  {
    j = unif_rand_int(0, i-1);
    if (j < k)
      samp[j] = i + INDEX_BASE;
  }
  
  PutRNGstate();
  
  return 0;
}



static inline void sort_insertion(const int len, len_t *const x)
{
  int i, j;
  len_t tmp;
  
  for (i=1; i<len; i++)
  {
    tmp = x[i];
    j = i - 1;
    
    while (j >= 0 && x[j] > tmp)
    {
      x[j+1] = x[j];
      j--;
    }
    
    x[j+1] = tmp;
  }
}


#endif
