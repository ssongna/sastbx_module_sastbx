#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/math/special_functions.hpp>


#include <sastbx/pr/regul.h>


namespace sastbx { namespace pr {
namespace boost_python{

  void wrap_regul_basic();

  namespace {

    void init_module()
    {
      using namespace boost::python;
      wrap_regul_basic();
    }

}}}}

BOOST_PYTHON_MODULE(sastbx_pr_ext)
{
  sastbx::pr::boost_python::init_module();
}
