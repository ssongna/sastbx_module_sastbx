#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <scitbx/math/zernike.h>
#include <sastbx/zernike_model/zernike_model.h>

namespace sastbx { namespace zmodel {
namespace boost_python{

  struct zmodel_wrapper
  {
    typedef zernike_model  < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("zernike_model", no_init)
        .def( init<
                   scitbx::math::zernike::nlm_array<double> const&,
                   scitbx::af::const_ref<double> const&,
                   double const& ,
                   int const&
                  >
             ((
                arg("Cnlm"),
                arg("q_array"),
                arg("r_max"),
                arg("n_max")
             ))
            )
        .def("calc_intensity", &w_t::calc_intensity)
        .def("calc_intensity_nnl", &w_t::calc_intensity_nnl)
        .def("calc_intensity_nlm", &w_t::calc_intensity_nlm)
      ;
    }
  };

  void
  wrap_zmodel()
  {
    zmodel_wrapper::wrap();
  }
//
//

}}}
