#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

#include <sastbx/zernike_model/zernike_model.h>

namespace sastbx { namespace zmodel {
namespace boost_python{

  void wrap_zmodel();
  namespace {

    void init_module()
    {
      using namespace boost::python;
      wrap_zmodel();
    }

}}}} // sastbx::zmodel::<anonymous>

BOOST_PYTHON_MODULE(sastbx_zmodel_ext)
{
  sastbx::zmodel::boost_python::init_module();
}
