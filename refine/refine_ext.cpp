#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/math/special_functions.hpp>

#include <sastbx/refine/elastic.h>
#include <sastbx/refine/elastic_rtb.h>

namespace sastbx { namespace refine {
namespace boost_python{

  void wrap_elastic();
  void wrap_elastic_rtb();
  namespace {

    void init_module()
    {
      using namespace boost::python;
      wrap_elastic();
      wrap_elastic_rtb();
    }

}}}} // sastbx::refine::<anonymous>

BOOST_PYTHON_MODULE(sastbx_refine_ext)
{
  sastbx::refine::boost_python::init_module();
}
