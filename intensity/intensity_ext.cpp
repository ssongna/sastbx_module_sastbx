#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/math/special_functions.hpp>

#include <sastbx/intensity/model.h>
#include <sastbx/intensity/scatlib.h>
#include <sastbx/intensity/debye_engine.h>
#include <sastbx/intensity/she_engine.h>
#include <sastbx/intensity/integrated.h>
#include <sastbx/intensity/besselR.h>

namespace sastbx { namespace intensity {
namespace boost_python{

  void wrap_model();
  void wrap_debye_engine();
  void wrap_she_engine();
  void wrap_integrated();
  void wrap_block_integrator();
  void wrap_besselr();
  void wrap_scattering_library();

  namespace {

    void init_module()
    {
      using namespace boost::python;
      wrap_model();
      wrap_debye_engine();
      wrap_she_engine();
      wrap_integrated();
      wrap_block_integrator();
      wrap_besselr();
      wrap_scattering_library();
    }

}}}} // sastbx::intensity::<anonymous>

BOOST_PYTHON_MODULE(sastbx_intensity_ext)
{
  sastbx::intensity::boost_python::init_module();
}
