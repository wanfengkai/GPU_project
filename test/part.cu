#include <stdio.h>
#include "part.h"

const int blocksize = 32;

__global__ void hello(int a)
{
    printf("Hello %d! My threadId is %d\n",a, threadIdx.x );
}

void launchkernel(int a)
{
   hello<<<1, blocksize>>>(a);
   cudaDeviceSynchronize();
}
