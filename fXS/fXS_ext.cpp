#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <sastbx/fXS/math_utils.h>
#include <sastbx/fXS/images.h>
#include <sastbx/fXS/ewald_sphere.h>
#include <sastbx/fXS/image_simulator_tools.h>
#include <sastbx/fXS/zernike_model_tools.h>
#include <sastbx/fXS/detector_tools.h>

namespace sastbx { namespace fXS {
namespace {

  struct image_wrapper
  {
    typedef image_cartesian  < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("image_cartesian", no_init)
        .def( init<
                   int const&,
                   int const&,
                   scitbx::af::const_ref< scitbx::vec3<double> >
                  >
             ((
                arg("np"),
                arg("n_max"),
                arg("image")
             ))
            )
        .def("calc_c2_array", &w_t::calc_c2_array)
        .def("get_c2_array", &w_t::get_c2_array)
        .def("get_normalized_sp", &w_t::get_normalized_sp)
      ;
    }
  };

// wrapper for Blq expansion coefficient calculation
// using zernike expansion method
  struct zernike_wrapper
  {
    typedef zernike_moment_variants< double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("zernike_moment_variants", no_init)
        .def( init<
                   scitbx::math::zernike::nlm_array<double> const&,
                   scitbx::af::const_ref< double >,
                   double const&,
                   int const&,
                   int const&
                  >
             ((
                arg("c_nlm"),
                arg("q_array"),
                arg("rmax"),
                arg("nmax"),
                arg("lmax")
             ))
            )
        .def("get_all_blq", &w_t::get_all_blq)
        .def("get_all_blq2", &w_t::get_all_blq2)
        .def("get_real_coef", &w_t::get_real_coef)
        .def("get_alm", &w_t::get_a_expansion)
        .def("get_ilm", &w_t::get_i_expansion)
      ;
    }
  };



  struct ewald_sphere_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<sastbx::fXS::ewald_sphere>("ewald_sphere", init<>() )
        .def(init<const double&, const double& >() )
        .def("set_wavelength",
             &sastbx::fXS::ewald_sphere::set_wavelength)
        .def("set_distance",
             &sastbx::fXS::ewald_sphere::set_distance)
        .def("get_distance",
             &sastbx::fXS::ewald_sphere::get_distance)
        .def("get_h",
             &sastbx::fXS::ewald_sphere::get_h)
        .def("h_to_q",
             &sastbx::fXS::ewald_sphere::h_to_q)
      ;
    }
  };

  struct image_composer_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<sastbx::fXS::image_composer>("image_composer", init<>() )
        .def(init<const scitbx::vec2<int>&,const scitbx::vec2<int>&,
             const double&,const sastbx::fXS::ewald_sphere& >() )
        .def("set_detector_size",
             &sastbx::fXS::image_composer::set_detector_size)
        .def("get_detector_size",
             &sastbx::fXS::image_composer::get_detector_size)
        .def("set_beam_center",
             &sastbx::fXS::image_composer::set_beam_center)
        .def("get_beam_center",
             &sastbx::fXS::image_composer::get_beam_center)
        .def("set_pixel_size",
             &sastbx::fXS::image_composer::set_pixel_size)
        .def("get_pixel_size",
             &sastbx::fXS::image_composer::get_pixel_size)
        .def("set_ewald_sphere",
             &sastbx::fXS::image_composer::set_ewald_sphere)
        .def("get_ewald_sphere",
             &sastbx::fXS::image_composer::get_ewald_sphere)
        .def("get_distance",
             &sastbx::fXS::image_composer::get_distance)
        .def("cache_h",
             &sastbx::fXS::image_composer::cache_h)
        .def("build_image",
             &sastbx::fXS::image_composer::build_image)
        .def("get_q",
             &sastbx::fXS::image_composer::get_q)
        .def("get_phi",
             &sastbx::fXS::image_composer::get_phi)
        ;
    }
  };

  struct c2_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<sastbx::fXS::c2>("c2", init<>() )
        .def("set_image_composer",&sastbx::fXS::c2::set_image_composer)
        .def("set_q_limits",&sastbx::fXS::c2::set_q_limits)
        .def("set_pixel_sizes",&sastbx::fXS::c2::set_pixel_sizes)
        .def("initialize",&sastbx::fXS::c2::initialize)
        .def("bin_mask",&sastbx::fXS::c2::bin_mask)
        .def("process_intensities",&sastbx::fXS::c2::process_intensities)
        .def("get_ring_q",&sastbx::fXS::c2::get_ring_q)
        .def("get_n_rings",&sastbx::fXS::c2::get_n_rings)
        .def("get_c2",&sastbx::fXS::c2::get_c2)
        ;
    }
  };

  struct detector_geometry_wrapper {
    static void wrap() {
      using namespace boost::python;
      class_<sastbx::fXS::detector_geometry>("detector_geometry",init<>() )
        .def("set_corner_position",
             &sastbx::fXS::detector_geometry::set_corner_position)
        .def("get_corner_position",
             &sastbx::fXS::detector_geometry::get_corner_position)
        .def("set_detector_size",
             &sastbx::fXS::detector_geometry::set_detector_size)
        .def("get_detector_size",
             &sastbx::fXS::detector_geometry::get_detector_size)
        .def("set_pixel_size",&sastbx::fXS::detector_geometry::set_pixel_size)
        .def("get_pixel_size",&sastbx::fXS::detector_geometry::get_pixel_size)
        ;
    }
  };

  struct image_base_wrapper {
    static void wrap() {
      using namespace boost::python;
      class_<sastbx::fXS::image_base>("image_base",init<>() )
        .def("set_ewald_sphere",&sastbx::fXS::image_base::set_ewald_sphere)
        .def("get_ewald_sphere",&sastbx::fXS::image_base::get_ewald_sphere)
        .def("set_detector_geometry",
             &sastbx::fXS::image_base::set_detector_geometry)
        .def("get_detector_geometry",
             &sastbx::fXS::image_base::get_detector_geometry)
        .def("reset",&sastbx::fXS::image_base::reset)
        .def("get_center_q",&sastbx::fXS::image_base::get_center_q)
        .def("get_center_phi",&sastbx::fXS::image_base::get_center_phi)
        .def("get_corner_h",&sastbx::fXS::image_base::get_corner_h)
        .def("integrate",&sastbx::fXS::image_base::integrate)
        ;
    }
  };

  struct c2_tiles_wrapper {
    static void wrap() {
      using namespace boost::python;
      class_<sastbx::fXS::c2_tiles>("c2_tiles",init<>() )
        .def("set_ewald_sphere",&sastbx::fXS::c2_tiles::set_ewald_sphere)
        .def("add_geometry",&sastbx::fXS::c2_tiles::add_geometry)
        .def("add_intensities",&sastbx::fXS::c2_tiles::add_intensities)
        .def("set_q_limits",&sastbx::fXS::c2_tiles::set_q_limits)
        .def("set_ring_pixel_sizes",
             &sastbx::fXS::c2_tiles::set_ring_pixel_sizes)
        .def("initialize",&sastbx::fXS::c2_tiles::initialize)
        .def("bin_mask",&sastbx::fXS::c2_tiles::bin_mask)
        .def("process_intensities",&sastbx::fXS::c2_tiles::process_intensities)
        .def("process_ring",&sastbx::fXS::c2_tiles::process_ring)
        .def("reset_intensities",&sastbx::fXS::c2_tiles::reset_intensities)
        .def("get_ring_q",&sastbx::fXS::c2_tiles::get_ring_q)
        .def("get_mean_ring_intensities",
             &sastbx::fXS::c2_tiles::get_mean_ring_intensities)
        .def("get_ring_intensities",&sastbx::fXS::c2_tiles::get_ring_intensities)
        .def("get_n_rings",&sastbx::fXS::c2_tiles::get_n_rings)
        .def("get_c2",&sastbx::fXS::c2_tiles::get_c2)
        .def("get_pixel_indices",&sastbx::fXS::c2_tiles::get_pixel_indices)
        ;
    }
  };

  struct fourier_legendre_series_wrapper {
    static void wrap() {
      using namespace boost::python;
      class_<sastbx::fXS::fourier_legendre_series>
        ("fourier_legendre_series",init<>() )
        .def("precompute_polynomials",
             &sastbx::fXS::fourier_legendre_series::precompute_polynomials)
        .def("write_polynomials",
             &sastbx::fXS::fourier_legendre_series::write_polynomials)
        .def("read_polynomials",
             &sastbx::fXS::fourier_legendre_series::read_polynomials)
        .def("P",&sastbx::fXS::fourier_legendre_series::P)
        .def("compute_coefficients",
             &sastbx::fXS::fourier_legendre_series::compute_coefficients)
        .def("compute_coefficient_variances",
             &sastbx::fXS::fourier_legendre_series::compute_coefficient_variances)
        .def("compute_coefficients_ww",
             &sastbx::fXS::fourier_legendre_series::compute_coefficients_ww)
        .def("compute_function",
             &sastbx::fXS::fourier_legendre_series::compute_function)
        ;
    }
  };

}

namespace boost_python {

  void wrap_images()
  {
    image_wrapper::wrap();
  }

  void wrap_zernike()
  {
    zernike_wrapper::wrap();
  }

  void wrap_ewald_sphere()
  {
    ewald_sphere_wrapper::wrap();
  }

  void wrap_image_composer()
  {
    image_composer_wrapper::wrap();
  }

  void wrap_c2()
  {
    c2_wrapper::wrap();
  }

  void wrap_functions()
  {
    boost::python::def("mean_and_variance_by_q",&mean_and_variance_by_q);
    boost::python::def("I_q",&I_q);
    boost::python::def("direct_sum_structure_factors",
                       &direct_sum_structure_factors);
    boost::python::def("solvent_image",&solvent_image);
    boost::python::def("nearest_neighbors",&nearest_neighbors);
    boost::python::def("solvent_accessible_area",&solvent_accessible_area);
    boost::python::def("bilinear_interpolation",&bilinear_interpolation);
    boost::python::def("bauer_spiral",&bauer_spiral);
    boost::python::def("multiple_poisson",&multiple_poisson);
    boost::python::def("set_negative_to_zero",&set_negative_to_zero);
    boost::python::def("apply_translation",&apply_translation);
    boost::python::def("get_binning",&get_binning);
    boost::python::def("prebin_image",&prebin_image);
    boost::python::def("bin_mask",&bin_mask);
    boost::python::def("bin_intensities",&bin_intensities);
    boost::python::def("autocorrelation",&autocorrelation);
  }

}}}

BOOST_PYTHON_MODULE(sastbx_fXS_ext)
{
  sastbx::fXS::boost_python::wrap_images();
  sastbx::fXS::boost_python::wrap_zernike();
  sastbx::fXS::boost_python::wrap_ewald_sphere();
  sastbx::fXS::boost_python::wrap_image_composer();
  sastbx::fXS::boost_python::wrap_c2();
  sastbx::fXS::detector_geometry_wrapper::wrap();
  sastbx::fXS::image_base_wrapper::wrap();
  sastbx::fXS::c2_tiles_wrapper::wrap();
  sastbx::fXS::fourier_legendre_series_wrapper::wrap();
  sastbx::fXS::boost_python::wrap_functions();
}
