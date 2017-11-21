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

#include <sastbx/intensity/model.h>
#include <sastbx/intensity/debye_engine.h>
#include <sastbx/intensity/she_engine.h>
#include <sastbx/intensity/scatlib.h>
#include <sastbx/intensity/besselR.h>
#include <sastbx/intensity/integrated.h>
#include <cctbx/xray/scatterer.h>

namespace sastbx { namespace intensity {

namespace{

  // scattering factor library
  struct scattering_library_wrapper
  {
    typedef scattering_library< double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("scattering_library", no_init)
        .def( init< scitbx::af::const_ref< double > const& > (( arg_("q_values") )) )
        .def( "q", &w_t::q )
        .def( "q_range", &w_t::q_range)
        .def( "load_scattering_info", &w_t::load_scattering_info)
        .def( "get_sfs", &w_t::get_sfs)
        .def( "get_sf_indx", &w_t::get_sf_indx)
        .def( "get_sf", &w_t::get_sf)
      ;
    }

  };
  //Integration Mapping
  struct integrated_wrapper
  {
    typedef integrated < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("integrated", no_init)
        .def( init<
                double const&,
                double const&,
                double const&,
                double const&,
                double const&,
                double const&,
                int const&
                >
        (( arg_("f_min"),
           arg_("f_max"),
           arg_("delta"),
           arg_("resl"),
           arg_("q_step"),
           arg_("q_max"),
           arg_("l")
        ))
        )
        .def( "Get_integral", &w_t::Get_integral )
        .def( "print_out", &w_t::print_out )
      ;
    }
  };


  struct block_integrator_wrapper
  {
    typedef block_integrator < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("block_integrator", no_init)
        .def( init<
                double const&,
                double const&,
                double const&,
                double const&
                >
        (( arg_("dmax"),
           arg_("width"),
           arg_("q_max"),
           arg_("q_step")
        ))
        )
        .def( "setup_arrays", &w_t::setup_arrays)
        .def( "get_intensity", &w_t::get_intensity)
      ;
    }
  };


  struct besselr_wrapper
  {
    typedef besselr < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("besselr", no_init)
        .def( init<
                double const&,
                double const&,
                int const&
                >
        (( arg_("q_step"),
           arg_("q_max"),
           arg_("l")
        ))
        )
        .def("print_out", &w_t::print_out)
      ;
    }

  };

  //Model
  struct model_wrapper
  {
    typedef model < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("model", no_init)
        .def( init<
                     scitbx::af::const_ref< scitbx::vec3<double> > const&,
                     scitbx::af::const_ref< double > const&,
                     scitbx::af::const_ref< double > const&,
                     scitbx::af::const_ref< double > const&,
                     scitbx::af::const_ref< std::string > const&,
                     sastbx::intensity::scattering_library< double > const&,
                     bool const&
                  >
          ((
             arg_("xyz"),
             arg_("radii"),
             arg_("b_values"),
             arg_("occupancy"),
             arg_("atom_types"),
             arg_("scat_lib"),
             arg_("B_factor_on")
          ))
        )
        .def("get_max_radius", &w_t::get_max_radius)
        .def("size", &w_t::size)
        .def("get_xyz", &w_t::get_xyz)
        ;
    }

  };

  // Debye engine
  struct debye_engine_wrapper
  {
    typedef debye_engine < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("debye_engine", no_init)
        .def( init<
                   sastbx::intensity::model <double> const&,
                   sastbx::intensity::scattering_library< double > const&
                  >
          ((
             arg_("model"),
             arg_("scat_lib")
          ))
        )
        .def("I", &w_t::I)
        ;
    }

  };


  // She engine
  struct she_engine_wrapper
  {
    typedef she_engine< double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("she_engine", no_init)
        .def( init<
                     sastbx::intensity::model <double> const&,
                     sastbx::intensity::scattering_library< double > const&,
                     int const&,
                     int const&,
                     double const&,
                     double const&,
                     double const&,
                     double const&,
                     double const&,
                     double const&
                  >
          ((
             arg_("model"),
             arg_("scat_lib_dummy"),
             arg_("max_i"),
             arg_("max_L"),
             arg_("f_resl"),
             arg_("q_step"),
             arg_("max_z"),
             arg_("delta"),
             arg_("rho"),
             arg_("d_rho")
          ))
        )
        .def( "I", &w_t::I )
        .def( "Iscale", &w_t::Iscale )
        .def( "Compute_Coefs", &w_t::Compute_Coefs )
        .def( "get_IA", &w_t::get_IA )
        .def( "get_IB", &w_t::get_IB )
        .def( "get_IC", &w_t::get_IC )
        .def( "Area_Volume", &w_t::Area_Volume )
        .def( "Area_Volume2", &w_t::Area_Volume2 )
        .def( "update_coord", &w_t::update_coordinate )
        .def( "update_solvent_params", &w_t::UpdateSolventParams )
	.def( "calc_spatial_correlation", &w_t::calc_spatial_correlation)
	.def( "get_spatial_correlation", &w_t::get_spatial_correlation)
	.def( "get_expansion_coef", &w_t::get_expansion_coef)
	.def( "get_all_coefs", &w_t::get_all_coefs)
	.def( "A_q", &w_t::A_q)
	.def( "I_q_omega_default", &w_t::I_q_omega_default)
	.def( "get_grid_xyz", &w_t::get_grid_xyz)
	.def( "I_q_omega", &w_t::I_q_omega)
	.def( "simulate_correlation", &w_t::simulate_correlation)
        ;
    }

  };



} // namespace <anonymous>

namespace boost_python{

  void wrap_model()
  {
   sastbx::intensity::model_wrapper::wrap();
  }

  void wrap_debye_engine()
  {
   sastbx::intensity::debye_engine_wrapper::wrap();
  }

  void wrap_she_engine()
  {
    sastbx::intensity::she_engine_wrapper::wrap();
  }

  void wrap_integrated()
  {
    sastbx::intensity::integrated_wrapper::wrap();
  }
  void wrap_block_integrator()
  {
    sastbx::intensity::block_integrator_wrapper::wrap();
  }

  void wrap_besselr()
  {
    sastbx::intensity::besselr_wrapper::wrap();
  }

  void wrap_scattering_library()
  {
    sastbx::intensity::scattering_library_wrapper::wrap();
  }


}

}}
