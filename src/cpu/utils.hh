#ifndef CLUSTRGAME_UTILS_H_
#define CLUSTRGAME_UTILS_H_


#define FALSE 0
#define TRUE 1

#define SQRT_EPS_FLT 1e-4f
#define SQRT_EPS_DBL 1e-8


static inline int all_equal(const int len, const float *const restrict x, const float *const restrict y)
{
  float mad = 0.0;
  
  for (int i=0; i<len; i++)
    mad += x[i] - y[i];
  
  mad /= ((float) len);
  
  return mad > SQRT_EPS_FLT ? FALSE : TRUE;
}

static inline int all_equal(const int len, const double *const restrict x, const double *const restrict y)
{
  double mad = 0.0;
  
  for (int i=0; i<len; i++)
    mad += x[i] - y[i];
  
  mad /= ((double) len);
  
  return mad > SQRT_EPS_DBL ? FALSE : TRUE;
}



static inline void add1(const int n, int *const x)
{
  for (int i=0; i<n; i++)
    x[i]++;
}


#endif
