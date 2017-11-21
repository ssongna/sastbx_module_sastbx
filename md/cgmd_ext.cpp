#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <sastbx/md/beads.h>

namespace sastbx { namespace cgmd {
namespace boost_python{

  void wrap_beads();
  namespace {

    void init_module()
    {
      using namespace boost::python;
      wrap_beads();
    }

}}}} // sastbx::cgmd::<anonymous>

BOOST_PYTHON_MODULE(sastbx_cgmd_ext)
{
  sastbx::cgmd::boost_python::init_module();
}
