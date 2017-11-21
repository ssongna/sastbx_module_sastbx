#ifndef DETECTOR_TOOLS_H

#include <math.h>
#include <vector>
#include <set>

#include <sastbx/fXS/ewald_sphere.h>
#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/math/basic_statistics.h>

namespace sastbx{
namespace fXS {

  // --------------------------------------------------------------------------
  class detector_geometry {

  public:
    detector_geometry();
    void set_corner_position(const scitbx::vec3<double>&);
    scitbx::vec3<double> get_corner_position();
    void set_detector_size(const scitbx::vec2<int>&);
    scitbx::vec2<int> get_detector_size();
    void set_pixel_size(const scitbx::vec2<double>&);
    scitbx::vec2<double> get_pixel_size();

  private:
    scitbx::vec3<double> corner_position;
    scitbx::vec2<int> detector_size;
    scitbx::vec2<double> pixel_size;
  };

  // --------------------------------------------------------------------------
  class image_base {

  public:
    image_base();
    void set_ewald_sphere(const ewald_sphere&);
    ewald_sphere get_ewald_sphere();
    void set_detector_geometry(const detector_geometry&);
    detector_geometry get_detector_geometry();
    void reset();
    scitbx::af::shared<double> get_center_q();
    scitbx::af::shared<double> get_center_phi();
    scitbx::af::shared<scitbx::vec3<double> > get_corner_h();
    scitbx::af::shared<double> integrate(const scitbx::af::const_ref<double>&);

  private:
    detector_geometry dg;
    ewald_sphere es;
    scitbx::af::shared<double> q;
    scitbx::af::shared<double> phi;
    scitbx::af::shared<scitbx::vec3<double> > h;
  };

  // --------------------------------------------------------------------------
  class c2_tiles {

  public:
    c2_tiles();
    void set_ewald_sphere(const ewald_sphere&);
    void add_geometry(const detector_geometry&);
    void add_intensities(const scitbx::af::const_ref<double>&);
    void set_q_limits(const double&, const double&);
    void set_ring_pixel_sizes(const int&, const int&);
    void initialize();
    scitbx::af::shared<int> bin_mask(const int&);
    void process_intensities();
    void process_ring(const int&);
    void reset_intensities();
    double get_ring_q(const int&);
    scitbx::af::shared<double> get_mean_ring_intensities();
    scitbx::af::shared<double> get_ring_intensities(const int&);
    int get_n_rings();
    scitbx::af::shared<double> get_c2(const int&);
    scitbx::af::shared<int> get_pixel_indices(const int&,const int&);

  private:
    std::vector<detector_geometry> tiles;
    std::vector<scitbx::af::const_ref<double> > intensities;
    scitbx::vec2<double> upper_left_xy, lower_right_xy;
    std::vector<scitbx::af::shared<double> > q;
    std::vector<scitbx::af::shared<double> > phi;
    image_base ib;
    double q_min, q_max;
    scitbx::af::shared<double> q0, phi0;
    int q_pixel_depth, phi_pixel_radius;
    std::vector<std::vector<std::vector<std::vector<int> > > >pixel_bins;
    std::vector<std::vector<std::vector<bool> > > use_bin;
    scitbx::af::shared<double> ring_q;
    scitbx::af::shared<double> mean_ring_intensities;
    std::vector<std::vector<double> > ring_intensities;
    std::vector<std::vector<double> > ring_c2;
    int get_i_from_q0(const double&);
  };

}
}
#endif // DETECTOR_TOOLS_H
