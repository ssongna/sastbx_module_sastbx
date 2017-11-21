#ifndef CUDA_FUNCTIONS_CUH
#define CUDA_FUNCTIONS_CUH

#include <sastbx/fXS/cuda_functions.h>
#include <cudatbx/cuda_base.cuh>

namespace sastbx {
namespace fXS {

  __global__ void structure_factor_kernel
    (const int*, const float*, const float*, const int,
     const float*, const int,
     const float*, const int,
     float*, float*);

  /* ==========================================================================
   */

}
}
#endif // CUDA_FUNCTIONS_CUH
