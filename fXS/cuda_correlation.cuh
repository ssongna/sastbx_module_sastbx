#ifndef CUDA_CORRELATION_CUH
#define CUDA_CORRELATION_CUH

#include <sastbx/fXS/cuda_correlation.h>
#include <cudatbx/cuda_base.cuh>
#include <cuComplex.h>

namespace sastbx {
namespace fXS {

  __global__ void add_images_kernel_streams
    (const cuDoubleComplex*, const int, const int, cuDoubleComplex*);
  scitbx::af::shared<std::complex<double> > cuda_add_images_streams
    (const scitbx::af::const_ref<std::complex<double> >&, const int&,
     const int&, const int&);

  __global__ void add_images_kernel
    (const cuDoubleComplex*, const int, const int, cuDoubleComplex*);

}
}
#endif // CUDA_CORRELATION_CUH
