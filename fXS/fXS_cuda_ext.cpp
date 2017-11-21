#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <sastbx/fXS/cuda_functions.h>
#include <sastbx/fXS/cuda_correlation.h>

namespace sastbx { namespace fXS {

namespace boost_python {

  struct cuda_direct_summation_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<sastbx::fXS::cuda_direct_summation>("cuda_direct_summation",init<>() )
        .def("add",&sastbx::fXS::cuda_direct_summation::add)
        .def("get_sum",&sastbx::fXS::cuda_direct_summation::get_sum)
        ;
    }
  };

  void wrap_functions()
  {
    boost::python::def("cuda_add_images_streams",&cuda_add_images_streams);
    boost::python::def("cuda_add_images",&cuda_add_images);
  }

}}

BOOST_PYTHON_MODULE(sastbx_fXS_cuda_ext)
{
  sastbx::fXS::boost_python::cuda_direct_summation_wrapper::wrap();
  sastbx::fXS::boost_python::wrap_functions();
}
}
