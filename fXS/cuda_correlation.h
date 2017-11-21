#ifndef CUDA_CORRELATION_H
#define CUDA_CORRELATION_H

#include <math.h>
#include <complex>

#include <scitbx/array_family/shared.h>

namespace sastbx {
namespace fXS {

  scitbx::af::shared<std::complex<double> > cuda_add_images_streams
    (const scitbx::af::const_ref<std::complex<double> >&, const int&,
     const int&, const int&);

  scitbx::af::shared<std::complex<double> > cuda_add_images
    (const scitbx::af::const_ref<std::complex<double> >&, const int&,
     const int&, const int&);

}
}
#endif // CUDA_CORRELATION_H
