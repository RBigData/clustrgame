#ifndef KMEANS_UTILS_H_
#define KMEANS_UTILS_H_


#define FALSE 0
#define TRUE 1

#define EPS 1e-8

static inline int all_equal(const int len, const double *const restrict x, const double *const restrict y)
{
  double mad = 0.0;
  
  for (int i=0; i<len; i++)
    mad += x[i] - y[i];
  
  mad /= ((double) len);
  
  return mad > EPS ? FALSE : TRUE;
}



static inline void add1(const int n, int *const x)
{
  for (int i=0; i<n; i++)
    x[i]++;
}


#endif
