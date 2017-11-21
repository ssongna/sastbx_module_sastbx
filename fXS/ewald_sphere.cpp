#include <sastbx/fXS/ewald_sphere.h>

namespace sastbx {
namespace fXS {

  /* ==========================================================================
     ewald_sphere Class

     This class implements an Ewald sphere for determining which reciprocal
     space vectors are visible.  The formulas are based on references from
     Rossman

       J. Appl. Cryst. (1979). 12, 225-238
       Acta Cryst. (1999). D55, 1631

     The origin of the reciprocal lattice is at (0,0,0) and the center of the
     Ewald sphere is at (0,0,-1/lambda).  The +z-axis is in the direction of
     the incident beam.

     Public Methods:
       ewald_sphere - 2 constructors
       set_wavelength - self-explanatory
       set_distance - self-explanatory
       get_distance - self-explanatory
       get_h - returns the reciprocal space vector for a detector coordinate
       h_to_q - converts h into q

     --------------------------------------------------------------------------
     Constructor

     Arguments:
       None
       lambda - wavelength (double)
       d      - distance from sample to detector (double)

     Returns:
       None

     Notes:
       Every method in the non-trivial constructor needs to be called before
         calling any other method
     --------------------------------------------------------------------------
  */
  sastbx::fXS::ewald_sphere::ewald_sphere() {}

  sastbx::fXS::ewald_sphere::ewald_sphere
  (const double& lambda,const double& d) {
    set_wavelength(lambda);
    set_distance(d);
  }

  void sastbx::fXS::ewald_sphere::set_wavelength(const double& lambda) {
    wavelength = lambda;
    sphere_center[0] = 0.0;
    sphere_center[1] = 0.0;
    sphere_center[2] = -1.0/wavelength;
  }

  void sastbx::fXS::ewald_sphere::set_distance(const double& d) {
    distance = d;
  }

  double sastbx::fXS::ewald_sphere::get_distance() {
    return distance;
  }

  /* --------------------------------------------------------------------------
     Converts a detector coordinate into a floating-point Miller index

     Arguments:
       xy - detector coordinates (scitbx::vec2<double>)

     Returns:
       Floating-point diffraction vector, h (scitbx::vec3<double>)

     Notes:
       Equation 3 from (1999) reference is used (method is based on the ratio
         of lengths for similar triangles)
     --------------------------------------------------------------------------
  */
  scitbx::vec3<double> sastbx::fXS::ewald_sphere::get_h
  (const scitbx::vec2<double>& xy) {
    double scale = wavelength * std::sqrt(xy[0] * xy[0] + xy[1] * xy[1] +
                                          distance * distance);
    scale = 1.0/scale;
    scitbx::vec3<double> h;
    h[0] = xy[0] * scale;
    h[1] = xy[1] * scale;
    h[2] = distance * scale;
    h += sphere_center;
    return h;
  }

  /* --------------------------------------------------------------------------
     Returns the q values for the pixels in the image

     Arguments:
       reciprocal space vector (scitbx::vec3<double>)

     Returns:
       q value (double)

     Notes:

               4 pi sin(theta)             2 sin(theta)
        |q| = -----------------     |h| = --------------
                   lambda                     lambda

     --------------------------------------------------------------------------
  */
  double sastbx::fXS::ewald_sphere::h_to_q(const scitbx::vec3<double>& h) {
    double stol = 0.5 * std::sqrt(h*h);
    double q = scitbx::constants::four_pi * stol;
    return q;
  }

  /* ==========================================================================
     end of ewald_sphere class
  */
}
}
