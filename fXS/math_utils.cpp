#include <sastbx/fXS/math_utils.h>

namespace sastbx {
namespace fXS {

  /* ==========================================================================
     Order of points is assumed to be f00, f10, f01, f11

       y1 -- f01 ----- f11
              |         |
              |         |
              |         |
       y0 -- f00 ----- f10
              |         |
              x0        x1

       (x,y,f)

     --------------------------------------------------------------------------
  */
  double bilinear_interpolation
  (const scitbx::vec2<double>& xy,
   const scitbx::af::const_ref<scitbx::vec3<double> >& f) {
    double x1x0 = f[1][0] - f[0][0];
    double y1y0 = f[2][1] - f[0][1];
    double x1x = f[1][0] - xy[0];
    double y1y = f[2][1] - xy[1];
    double xx0 = xy[0] - f[0][0];
    double yy0 = xy[1] - f[0][1];
    double f00 = f[0][2];
    double f10 = f[1][2];
    double f01 = f[2][2];
    double f11 = f[3][2];
    double f_i = (f00*x1x*y1y + f10*xx0*y1y + f01*x1x*yy0 + f11*xx0*yy0)/\
      (x1x0*y1y0);
    return f_i;
  }

  /* ==========================================================================
     Constructs a set of points uniformly distributed along a spherical spiral,
     'uniform' meaning the points are equidistant from one another (area
     represented by each point is the same)

     JOURNAL OF GUIDANCE, CONTROL, AND DYNAMICS
     Vol. 23, No. 1, January-February 2000

     "Distribution of Points on a Sphere with Application to Star Catalogs"
     Robert Bauer

     Algorithm 4

     --------------------------------------------------------------------------
  */
  scitbx::af::shared<scitbx::vec3<double> > bauer_spiral(const int& n) {

    scitbx::af::shared<scitbx::vec3<double> > points(n);
    double L = std::sqrt(scitbx::constants::pi * n);
    double x, y, z, theta, phi;
    for (int k=1; k<(n+1); k++) {
      z = 1.0 - (2.0*k - 1.0)/n;
      phi = std::acos(z);
      theta = L * phi;
      x = std::sin(phi) * std::cos(theta);
      y = std::sin(phi) * std::sin(theta);
      points[k-1][0] = x;
      points[k-1][1] = y;
      points[k-1][2] = z;
    }
    return points;
  }

  /* ==========================================================================
     wrapper for boost::random::poisson_distribution
     --------------------------------------------------------------------------
  */
  scitbx::af::shared<double> multiple_poisson
  (const scitbx::af::const_ref<double>& l,const long& seed) {
    double cutoff = 700.0;
    boost::random::mt19937 random_number(static_cast<unsigned>(seed));
    scitbx::af::shared<double> result(l.size());
    for (int i=0; i<l.size(); i++) {
      if (l[i] <= cutoff) {
        /*
        boost::random::poisson_distribution<long,double> distribution(l[i]);
        boost::variate_generator<boost::random::mt19937&,
                                 boost::random::poisson_distribution<long,double> >
          p(random_number,distribution);
        result[i] = p();
        */
        boost::random::uniform_01<double> distribution;
        boost::variate_generator<boost::random::mt19937&,
                                 boost::random::uniform_01<double> >
          u(random_number,distribution);
        double L = std::exp(-l[i]);
        long k = 0;
        double p = 1.0;
        while (p > L) {
          k++;
          p = p * u();
        }
        result[i] = double(k - 1);
      }
      else {
        boost::random::normal_distribution<double>
          distribution(l[i],std::sqrt(std::fabs(l[i])));
        boost::variate_generator<boost::random::mt19937&,
                                 boost::random::normal_distribution<double> >
          p(random_number,distribution);
        result[i] = p();
        if (result[i] < 0.0) {
          result[i] = 0.0;
        }
      }
    }
    return result;
  }

  /* ==========================================================================
     sets any negative elements in an array to zero
     --------------------------------------------------------------------------
  */
  scitbx::af::shared<double> set_negative_to_zero
  (const scitbx::af::const_ref<double>& a) {
    scitbx::af::shared<double> b(a.size());
    for (int i=0; i<a.size(); i++) {
      if (a[i] < 0.0) {
        b[i] = 0.0;
      }
      else {
        b[i] = a[i];
      }
    }
    return b;
  }

  /* ==========================================================================
     Fourier-Legendre Series
     --------------------------------------------------------------------------
  */
  sastbx::fXS::fourier_legendre_series::fourier_legendre_series() {
    ready = false;
    highest_order_n = -1;
    polynomial_size = -1;
  }

  void sastbx::fXS::fourier_legendre_series::precompute_polynomials
  (const int& n, const int& s) {

    // reset storage
    highest_order_n = n+1;
    polynomial_size = s;
    polynomials.clear();
    polynomials.resize(highest_order_n * polynomial_size);
    dydx.clear();
    dydx.resize(highest_order_n*(polynomial_size - 1));

    // calculate polynomials and slopes
    dx = 2.0/(polynomial_size - 1);
    x.clear();
    x.resize(polynomial_size);
    x[0] = -1.0;
    for (int i=1; i<polynomial_size; i++) {
      x[i] = x[i-1] + dx;
    }
    x[polynomial_size - 1] = 1.0;
    int p_shift, d_shift;
    for (int order=0; order<highest_order_n; order++) {
      p_shift = order*polynomial_size;
      for (int i=0; i<polynomial_size; i++) {
        polynomials[p_shift + i] = boost::math::legendre_p(order,x[i]);
      }
      d_shift = order*(polynomial_size - 1);
      for (int i=0; i<(polynomial_size - 1); i++) {
        dydx[d_shift + i] = (polynomials[p_shift + i + 1] -
                             polynomials[p_shift + i])/(x[i+1] - x[i]);
      }
    }

    ready = true;
  }

  void sastbx::fXS::fourier_legendre_series::write_polynomials
  (const std::string& filename) {
    SCITBX_ASSERT( ready );

    FILE * f = fopen(filename.c_str(),"wb");

    // write dimensions
    fwrite(&highest_order_n, 1, sizeof(int), f);
    fwrite(&polynomial_size, 1, sizeof(int), f);

    // write data
    fwrite(&(polynomials[0]), polynomials.size(), sizeof(double), f);
    fwrite(&(dydx[0]), dydx.size(), sizeof(double), f);
    fwrite(&(x[0]), x.size(), sizeof(double), f);

    fclose(f);
  }

  void sastbx::fXS::fourier_legendre_series::read_polynomials
  (const std::string& filename) {
    FILE * f = fopen(filename.c_str(),"rb");

    // read dimensions
    fread(&highest_order_n, 1, sizeof(int), f);
    fread(&polynomial_size, 1, sizeof(int), f);

    // reset storage and read data
    polynomials.clear();
    polynomials.resize(highest_order_n * polynomial_size);
    fread(&(polynomials[0]), polynomials.size(), sizeof(double), f);

    dydx.clear();
    dydx.resize(highest_order_n * (polynomial_size - 1));
    fread(&(dydx[0]), dydx.size(), sizeof(double), f);

    x.clear();
    x.resize(polynomial_size);
    fread(&(x[0]), x.size(), sizeof(double), f);

    fclose(f);

    dx = 2.0/(polynomial_size - 1);

    ready = true;
  }

  double sastbx::fXS::fourier_legendre_series::P
  (const int& order, const double& x_i) {
    // check bounds
    if (x_i < -1.0) {
      return polynomials[order*polynomial_size];
    }
    if (x_i > 1.0) {
      return polynomials[order*polynomial_size + polynomial_size - 1];
    }

    // find bounds for linear interpolation
    int i0 = int(floor((x_i + 1) / dx));
    if (i0 == (polynomial_size - 1)) {
      i0 = polynomial_size - 2;
    }
    int i1 = i0 + 1;
    double x0 = x[i0];
    double y0 = polynomials[order*polynomial_size + i0];

    return (y0 + (x_i - x0) * dydx[order*(polynomial_size - 1) + i0]);
  }

  scitbx::af::shared<double>
  sastbx::fXS::fourier_legendre_series::compute_coefficients
  (const int& order, const scitbx::af::const_ref<double>& f,
   const scitbx::af::const_ref<double>& x_i) {
    SCITBX_ASSERT( ready );
    SCITBX_ASSERT( f.size() == x_i.size() );
    SCITBX_ASSERT( order < highest_order_n );
    scitbx::af::shared<double> result(order+1,0.0);
    for (int i=0; i<order+1; i++) {
      // add contributions from the ends
      if (x_i[0] > -1.0) {
        result[i] += f[0] * P(i,x_i[0]) * (x_i[0] + 1.0);
      }
      if (x_i[-1] < 1.0) {
        result[i] += f[-1] * P(i,x_i[-1]) * (1.0 - x_i[-1]);
      }
      // use trapezoidal rule for integrating everything in the middle
      for (int j=0; j<x_i.size()-1; j++) {
        result[i] += ( 0.5 * ( f[j+1] * P(i,x_i[j+1]) + f[j] * P(i,x_i[j]) ) *
                       (x_i[j+1] - x_i[j]) );
      }
      result[i] = (i + 0.5) * result[i];
    }
    return result;
  }

  scitbx::af::shared<double>
  sastbx::fXS::fourier_legendre_series::compute_coefficient_variances
  (const int& order, const scitbx::af::const_ref<double>& f,
   const scitbx::af::const_ref<double>& x_i) {
    SCITBX_ASSERT( ready );
    SCITBX_ASSERT( f.size() == x_i.size() );
    SCITBX_ASSERT( order < highest_order_n );
    scitbx::af::shared<double> result(order+1,0.0);
    double dx_i,p;
    for (int i=0; i<order+1; i++) {
      p = P(i,x_i[0]);
      dx_i = 1.5 * (x_i[1] - x_i[0]);
      result[i] = f[0] * p * p * dx_i * dx_i;
      p = P(i,x_i[-1]);
      dx_i = 1.5 * (x_i[-1] - x_i[-2]);
      result[i] += f[-1] * p * p * dx_i * dx_i ;
      for (int j=1; j<x_i.size()-1; j++) {
        p = P(i,x_i[j]);
        dx_i = 0.5 * (x_i[j+1] - x_i[j-1]);
        result[i] += f[j] * p * p * dx_i * dx_i;
      }
      p = i + 0.5;
      result[i] = p * p * result[i];
    }
    return result;
  }

  scitbx::af::shared<double>
  sastbx::fXS::fourier_legendre_series::compute_coefficients_ww
  (const int& order, const scitbx::af::const_ref<double>& f,
   const scitbx::af::const_ref<double>& x_i,
   const scitbx::af::const_ref<double>& w_i ) {
    SCITBX_ASSERT( ready );
    SCITBX_ASSERT( f.size() == x_i.size() );
    SCITBX_ASSERT( order < highest_order_n );
    scitbx::af::shared<double> result(order+1,0.0);
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<x_i.size(); j++) {
        result[i] += f[j] * P(i,x_i[j]) * w_i[j];
      }
      result[i] = (i + 0.5) * result[i];
    }
    return result;
  }

  scitbx::af::shared<double>
  sastbx::fXS::fourier_legendre_series::compute_function
  (const scitbx::af::const_ref<double>& coefficients,
   const scitbx::af::const_ref<double>& x_i) {
    SCITBX_ASSERT( ready );
    SCITBX_ASSERT( coefficients.size() <= highest_order_n );
    scitbx::af::shared<double> result(x_i.size(),0.0);
    for (int i=0; i<x_i.size(); i++) {
      for (int j=0; j<coefficients.size(); j++) {
        result[i] += coefficients[j] * P(j,x_i[i]);
      }
    }
    return result;
  }

  /* --------------------------------------------------------------------------
   */

}
}
