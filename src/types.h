#ifndef KMEANS_TYPES_H_
#define KMEANS_TYPES_H_


#ifdef _cplusplus
#include <cstdint>
#else
#include <stdint.h>
#endif

#define MPI_LEN_T MPI_UINT64_T
#define MPI_LEN_LOCAL_T MPI_INT

typedef uint64_t len_t;
typedef int len_local_t;


#endif
