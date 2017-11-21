#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

#include <sastbx/rigid_body/wigner3j.h>
#include <sastbx/rigid_body/rigid_body.h>

namespace sastbx { namespace rigid_body {
namespace boost_python{

  void wrap_wigner3j_fast();
  void wrap_wigner3j();
  void wrap_wigner3j_zero();
  void wrap_rigidbody();
  namespace {

    void init_module()
    {
      using namespace boost::python;
      wrap_wigner3j_fast();
      wrap_wigner3j();
      wrap_wigner3j_zero();
      wrap_rigidbody();
    }

}

}}} // sastbx::rigid_body::<anonymous>

BOOST_PYTHON_MODULE(sastbx_rb_ext)
{
  sastbx::rigid_body::boost_python::init_module();
}
