#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/math/special_functions.hpp>

#include <sastbx/data_reduction/data_reduction.h>

namespace sastbx { namespace data_reduction {
namespace boost_python{

  void wrap_perpendicular_detector();
  void wrap_image_mask();
  void wrap_radial_integrator();
namespace {

  void init_module()
  {
    using namespace boost::python;
    wrap_perpendicular_detector();
    wrap_image_mask();
    wrap_radial_integrator();
  }

}}}} // sastbx::data_reduction::<anonymous>

BOOST_PYTHON_MODULE(sastbx_data_reduction_ext)
{
  sastbx::data_reduction::boost_python::init_module();
}
