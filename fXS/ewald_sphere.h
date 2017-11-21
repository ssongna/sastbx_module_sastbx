#ifndef EWALD_SPHERE_H
#define EWALD_SPHERE_H

#include <math.h>

#include <cctbx/miller.h>
#include <cctbx/uctbx.h>
#include <scitbx/constants.h>
#include <scitbx/mat3.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace sastbx {
namespace fXS {

  class ewald_sphere {

  public:
    ewald_sphere();
    ewald_sphere(const double&, const double&);
    void set_wavelength(const double&);
    void set_distance(const double&);
    double get_distance();
    scitbx::vec3<double> get_h(const scitbx::vec2<double>&);
    double h_to_q(const scitbx::vec3<double>&);

  private:
    double wavelength;
    scitbx::vec3<double> sphere_center;
    double distance;
  };

}
}
#endif // EWALD_SPHERE_H
