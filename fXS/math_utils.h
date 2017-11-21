#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <math.h>
#include <boost/boost/math/special_functions/legendre.hpp>
#include <boost/boost/random/mersenne_twister.hpp>
#include <boost/boost/random/normal_distribution.hpp>
#include <boost/boost/random/poisson_distribution.hpp>
#include <boost/boost/random/uniform_01.hpp>
#include <boost/boost/random/variate_generator.hpp>
#include <string>

#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>

namespace sastbx {
namespace fXS {

  // --------------------------------------------------------------------------

  double bilinear_interpolation
    (const scitbx::vec2<double>&,
     const scitbx::af::const_ref<scitbx::vec3<double> >&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<scitbx::vec3<double> > bauer_spiral(const int&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<double> multiple_poisson
    (const scitbx::af::const_ref<double>&,const long&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<double> set_negative_to_zero
    (const scitbx::af::const_ref<double>&);

  // --------------------------------------------------------------------------

  class fourier_legendre_series {

  public:
    fourier_legendre_series();
    void precompute_polynomials(const int&, const int&);
    void write_polynomials(const std::string&);
    void read_polynomials(const std::string&);
    double P(const int&, const double&);
    scitbx::af::shared<double> compute_coefficients
      (const int&, const scitbx::af::const_ref<double>&,
       const scitbx::af::const_ref<double>&);
    scitbx::af::shared<double> compute_coefficient_variances
      (const int&, const scitbx::af::const_ref<double>&,
       const scitbx::af::const_ref<double>&);
    scitbx::af::shared<double> compute_coefficients_ww
      (const int&, const scitbx::af::const_ref<double>&,
       const scitbx::af::const_ref<double>&,
       const scitbx::af::const_ref<double>&);
    scitbx::af::shared<double> compute_function
      (const scitbx::af::const_ref<double>&,
       const scitbx::af::const_ref<double>&);

  private:
    bool ready;
    int highest_order_n;
    int polynomial_size;
    scitbx::af::shared<double> polynomials;
    scitbx::af::shared<double> x;
    scitbx::af::shared<double> dydx;
    double dx;

  };

  // --------------------------------------------------------------------------
}
}
#endif // MATH_UTILS_H
