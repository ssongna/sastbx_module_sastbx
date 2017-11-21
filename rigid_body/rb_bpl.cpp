#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/histogram.h>

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

#include <sastbx/rigid_body/wigner3j.h>
#include <sastbx/rigid_body/rigid_body.h>

namespace sastbx { namespace rigid_body {

namespace{
  //fast wigner3j
  struct wigner3j_fast_wrapper
  {
    typedef wigner3j_fast < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("wigner3j_fast", no_init)
        .def( init< int const&
                  >
             (( arg_("max")
             ))
            )
        .def("compute", &w_t::compute)
      ;
    }

  };

  // general wigner3j
  struct wigner3j_wrapper
  {
    typedef wigner3j < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("wigner3j", no_init)
        .def( init< int const&,
                    int const&,
                    int const&,
                    int const&,
                    int const&,
                    int const&
                  >
             (( arg_("j1"),
                arg_("j2"),
                arg_("j3"),
                arg_("m1"),
                arg_("m2"),
                arg_("m3")
             ))
            )
        .def("check", &w_t::check)
        .def("get_value", &w_t::get_value)
      ;
    }

  };

// with all m's = 0
  struct wigner3j_zero_wrapper
  {
    typedef wigner3j_zero < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("wigner3j_zero", no_init)
        .def( init< int const&,
                    int const&,
                    int const&
                  >
             (( arg_("j1"),
                arg_("j2"),
                arg_("j3")
             ))
            )
        .def("get_value", &w_t::get_value)
      ;
    }

  };

// rigidbody class
  struct rigidbody_wrapper
  {
    typedef rigidbody < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("rigidbody", no_init)
        .def( init< scitbx::af::const_ref< scitbx::vec3 <double > > const&,
		    scitbx::af::const_ref< int >    const&,
		    double const&, 
		    int const&,
		    bool const&
		  >
	      (( arg_("coords"),
		 arg_("sample_indices"),
		 arg_("dmax"),
		 arg_("max_i"),
		 arg_("build_grid")
	      ))
	     )
	.def("reset",&w_t::reset)
	.def("get_hist",&w_t::get_histogram)
	.def("that_hist",&w_t::that_histogram)
	.def("rotate_translate",&w_t::rotate_translate)
	.def("rotate_only",&w_t::rotate_only)
	.def("rotate_around_one_point",&w_t::rotate_around_one_point)
	.def("rotate_around_two_point",&w_t::rotate_around_two_point)
	.def("translate_after_rotation",&w_t::translate_only)
	.def("rotate_translate",&w_t::rotate_translate_q)
	.def("rotate_translatev",&w_t::rotate_translate_v)
	.def("get_crd",&w_t::get_crd)
	.def("get_clash",&w_t::get_clash)
	.def("get_hist_var",&w_t::get_hist_var)
	.def("get_grid", &w_t::get_grid)
	.def("range_f_w", &w_t::range_f_w)
	.def("center", &w_t::get_center)
	.def("get_surface_atoms", &w_t::get_surface_atoms)
	.def("set_weights", &w_t::set_weights)
	.def("update_contact_list", &w_t::update_contact_list)
	.def("calc_contact_hist", &w_t::calc_contact_hist)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python{

  void wrap_wigner3j_fast()
  {
   sastbx::rigid_body::wigner3j_fast_wrapper::wrap();
  }

  void wrap_wigner3j()
  {
   sastbx::rigid_body::wigner3j_wrapper::wrap();
  }

  void wrap_wigner3j_zero()
  {
   sastbx::rigid_body::wigner3j_zero_wrapper::wrap();
  }

  void wrap_rigidbody()
  {
   sastbx::rigid_body::rigidbody_wrapper::wrap();
  }

}

}}
