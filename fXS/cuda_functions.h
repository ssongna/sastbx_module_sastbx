#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H

// std::system includes
#include <cmath>
#include <ctime>

#include <complex>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <math.h>
#include <boost/optional.hpp>
#include <cctbx/eltbx/xray_scattering.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>

namespace sastbx {
namespace fXS {

  class cuda_direct_summation {

  public:
    cuda_direct_summation();
    ~cuda_direct_summation();
    void add(const scitbx::af::const_ref<std::string>&,
             const scitbx::af::const_ref<scitbx::vec3<double> >&,
             const scitbx::af::const_ref<double>&,
             const scitbx::af::const_ref<scitbx::vec3<double> >&,
             const scitbx::af::const_ref<double>&,
             const scitbx::af::const_ref<scitbx::vec3<double> >&,
             const cctbx::xray::scattering_type_registry&);
    scitbx::af::shared<std::complex<double> > get_sum();

  private:
    int sf_size;
    float* sf_real, * sf_imag;
  };

  /* ==========================================================================
   */

}
}
#endif // CUDA_FUNCTIONS_H
