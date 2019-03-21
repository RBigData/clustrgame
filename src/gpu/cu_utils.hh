#ifndef CLUSTRGAME_CU_UTILS_H_
#define CLUSTRGAME_CU_UTILS_H_


#include <cuda_runtime.h>

#define TPB 512

#define CUFREE(x) {if(x)cudaFree(x);}
#define PRINT_CUDA_ERROR() printf("%s\n", cudaGetErrorString(cudaGetLastError()));


#endif
