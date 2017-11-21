#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/boost_python/is_polymorphic_workaround.h>


#include <sastbx/data_reduction/data_reduction.h>
#include <cctbx/xray/scatterer.h>

namespace sastbx { namespace data_reduction {

namespace{

  // simple detector
  struct perpendicular_detector_wrapper
  {
    typedef perpendicular_detector< double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("perpendicular_detector", no_init)
        .def(init< int const&, int const&, double const&, double const&, double const&, double const&, double const&, double const& >
             ((arg_("nx"),
               arg_("ny"),
               arg_("scalex"),
               arg_("scaley"),
               arg_("orgx"),
               arg_("orgy"),
               arg_("distance"),
               arg_("wavelength")
               )))
        .enable_pickling()
        .def( "compute_bin_indices", &w_t::compute_bin_indices )
        .def( "q_low_bin_values", &w_t::q_low_bin_values)
        .def( "q_mean_bin_values", &w_t::q_mean_bin_values)
        .def( "q_range", &w_t::q_range)
        .def( "nx", &w_t::nx)
        .def( "ny", &w_t::ny)
        .def( "scale_x", &w_t::scale_x)
        .def( "scale_y", &w_t::scale_y)
        .def( "origin_x", &w_t::origin_x)
        .def( "origin_y", &w_t::origin_y)
        .def( "distance", &w_t::distance)
        .def( "wavelength", &w_t::wavelength)
        .def( "q", &w_t::q)
        .def( "binning_setup", &w_t::binning_setup)
        ;
    }

  };

  // beamstop, edge and shadow masks
  struct image_mask_wrapper
  {
    typedef image_mask< double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("image_mask", no_init)
        .def(init< int const&,
                   int const& >
             (( arg_("nx"),
                arg_("ny")
               )))
         .enable_pickling()
         .def("mask", &w_t::mask)
         .def("add", &w_t::add)
         .def("add_circle", &w_t::add_circle)
         .def("add_polygon", &w_t::add_polygon)
         .def("detect_dead_pixels", &w_t::detect_dead_pixels)
         .def("mask_invert", &w_t::mask_invert)
        ;
    }

  };


  // an integrator
  struct radial_integrator_wrapper
  {
    typedef radial_integrator< double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("radial_integrator", no_init)
        .def(init< perpendicular_detector<double> const&,
                   scitbx::af::const_ref<int> const&,
                   scitbx::af::const_ref<int> const& >
                 (( arg_("detector"),
                    arg_("mask"),
                    arg_("shadow")
                 )))
         .enable_pickling()
         .def("integrate", &w_t::integrate)
         .def("mean_intensity", &w_t::mean_intensity)
         .def("variance_intensity", &w_t::variance_intensity)
         ;
    }

  };




} // namespace <anonymous>

namespace boost_python{

  void wrap_perpendicular_detector()
  {
    sastbx::data_reduction::perpendicular_detector_wrapper::wrap();
  }

  void wrap_image_mask()
  {
    sastbx::data_reduction::image_mask_wrapper::wrap();
  }

  void wrap_radial_integrator()
  {
    sastbx::data_reduction::radial_integrator_wrapper::wrap();
  }




}

}}
