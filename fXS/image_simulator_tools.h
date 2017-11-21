#ifndef IMAGE_SIMULATOR_TOOLS_H
#define IMAGE_SIMULATOR_TOOLS_H

#include <algorithm>
#include <complex>
#include <map>
#include <math.h>
#include <vector>

#include <cctbx/eltbx/xray_scattering.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <scitbx/constants.h>
#include <scitbx/math/utils.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>
#include <sastbx/fXS/math_utils.h>
#include <sastbx/fXS/ewald_sphere.h>

namespace sastbx {
namespace fXS {

  // --------------------------------------------------------------------------

  class image_composer {

  public:
    image_composer();
    image_composer(const scitbx::vec2<int>&,const scitbx::vec2<int>&,
                   const double&,const sastbx::fXS::ewald_sphere&);
    void set_detector_size(const scitbx::vec2<int>&);
    scitbx::vec2<int> get_detector_size();
    void set_beam_center(const scitbx::vec2<int>&);
    scitbx::vec2<int> get_beam_center();
    void set_pixel_size(const double&);
    double get_pixel_size();
    void set_ewald_sphere(const sastbx::fXS::ewald_sphere&);
    sastbx::fXS::ewald_sphere get_ewald_sphere();
    double get_distance();
    scitbx::af::shared<scitbx::vec3<double> > cache_h();
    scitbx::af::shared<double> build_image(const scitbx::af::const_ref<double>&);
    scitbx::af::shared<double> get_q();
    scitbx::af::shared<double> get_phi();

  private:
    scitbx::vec2<int> detector_size;
    scitbx::vec2<int> beam_center;
    double pixel_size;
    scitbx::af::shared<scitbx::vec2<int> > pixel_xy;
    sastbx::fXS::ewald_sphere ewald_sphere;
    scitbx::af::shared<scitbx::vec3<double> > h_cache;
  };

  // --------------------------------------------------------------------------

  scitbx::af::shared<scitbx::vec2<double> > mean_and_variance_by_q
    (const double&,const double&,const double&,const int&,
     const scitbx::af::const_ref<double>&,
     const scitbx::af::const_ref<double>&);

  // --------------------------------------------------------------------------

  double I_q(const scitbx::af::const_ref<double>&,
             const scitbx::af::const_ref<double>&,
             const double&,const double&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<std::complex<double> > direct_sum_structure_factors
    (const scitbx::af::const_ref<std::string>&,
     const scitbx::af::const_ref<scitbx::vec3<double> >&,
     const scitbx::af::const_ref<double>&,
     const scitbx::af::const_ref<scitbx::vec3<double> >&,
     const cctbx::xray::scattering_type_registry&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<double> solvent_image
    (const scitbx::af::const_ref<double>&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<int> nearest_neighbors
    (const scitbx::af::const_ref<scitbx::vec3<double> >&,
     const scitbx::af::const_ref<double>&, const int&, const double&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<int> solvent_accessible_area
    (const scitbx::af::const_ref<scitbx::vec3<double> >&,
     const scitbx::af::const_ref<double>&,
     const scitbx::af::const_ref<int>&, const double&, const int&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<std::complex<double> > apply_translation
    (const scitbx::af::const_ref<std::complex<double> >&,
     const scitbx::af::const_ref<scitbx::vec3<double> >&,
     const scitbx::vec3<double>&);


  // --------------------------------------------------------------------------
  class c2 {

  public:
    c2();
    void set_image_composer(const sastbx::fXS::image_composer&);
    void set_q_limits(const double, const double);
    void set_pixel_sizes(const int, const int);
    void initialize(const bool);
    scitbx::af::shared<int> bin_mask();
    void process_intensities(const scitbx::af::const_ref<double>&);
    double get_ring_q(const int&);
    int get_n_rings();
    scitbx::af::shared<double> get_c2(const int&);

  private:
    sastbx::fXS::image_composer ic;
    scitbx::vec2<int> beam_center;
    scitbx::vec2<int> detector_size;
    scitbx::af::shared<double> q, phi;
    double q_min, q_max;
    int q_pixel_depth, phi_pixel_radius;
    int get_q_index(const double&);
    std::vector<std::vector<std::vector<int> > > pixel_bins;
    scitbx::af::shared<double> ring_q;
    std::vector<std::vector<double> > ring_c2;

  };

  // --------------------------------------------------------------------------

  scitbx::vec3<double> get_binning
    (const scitbx::vec2<int>&, const scitbx::vec2<int>&,
     const scitbx::af::const_ref<double>&,
     const scitbx::af::const_ref<double>&,
     const double&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<int> prebin_image
    (const scitbx::vec3<double>&,
     const scitbx::af::const_ref<double>&,
     const scitbx::af::const_ref<double>&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<int> bin_mask
    (const scitbx::vec3<double>&,
     const scitbx::af::const_ref<double>&,
     const scitbx::af::const_ref<double>&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<double> bin_intensities
    (const scitbx::vec3<double>&,
     const scitbx::af::const_ref<double>&,
     const scitbx::af::const_ref<double>&,
     const scitbx::af::const_ref<double>&);

  // --------------------------------------------------------------------------

  scitbx::af::shared<double> autocorrelation
    (const scitbx::af::const_ref<double>&);

  // --------------------------------------------------------------------------
}
}
#endif // IMAGE_SIMULATOR_TOOLS_H
