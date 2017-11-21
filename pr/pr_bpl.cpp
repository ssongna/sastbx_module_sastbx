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

#include <sastbx/intensity/integrated.h>
#include <sastbx/pr/regul.h>

namespace sastbx { namespace pr {
  namespace{


    //Integration Mapping
    struct regul_basic_wrapper
    {
      typedef regul_basic < double > w_t;
      static void
      wrap()
      {
        using namespace boost::python;
        class_<w_t>("regul_basic", no_init)
          .def( init<
                  double const&,
                  int const&,
                  double const&,
                  double const&
                  >
          (( arg_("d_max"),
             arg_("n_points"),
             arg_("q_max"),
             arg_("q_step")
          ))
         )
          .def( "load_data", &w_t::load_data )
          .def( "load_p", &w_t::load_p )
          .def( "load_scales", &w_t::load_scales )
          .def( "r", &w_t::r )
          .def( "get_intensity", &w_t::get_intensity )
          .def( "get_delta", &w_t::get_delta )
          .def( "dt_dp", &w_t::dt_dp )
          .def( "hessian_t", &w_t::hessian_t )
          .def( "dt_ds", &w_t::dt_ds )
          .def( "ddt_dds", &w_t::ddt_dds )
          .def( "t_int", &w_t::t_int )
          .def( "t_regul", &w_t::t_regul )
          .def( "dregul_dp", &w_t::dregul_dp )
          .def( "hessian_regul", &w_t::hessian_regul)
        ;
      }
    };


  } // namespace <anonymous>

  namespace boost_python{
     void wrap_regul_basic()
     {
     sastbx::pr::regul_basic_wrapper::wrap();
     }
  }
}}
