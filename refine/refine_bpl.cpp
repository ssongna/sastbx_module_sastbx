#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/histogram.h>

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

#include <sastbx/refine/elastic.h>
#include <sastbx/refine/elastic_rtb.h>

namespace sastbx { namespace refine {

namespace{

  // CA based ENM //
  struct elastic_wrapper
  {
    typedef elastic < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("elastic", no_init)
        .def( init< scitbx::af::const_ref< scitbx::vec3 <double > > const&,
                    scitbx::af::const_ref< int >    const&,
                    double const&,
                    double const&
                  >
             (( arg_("xyz"),
                arg_("ca_indx"),
                arg_("cutoff"),
                arg_("scale_factor")
             ))
            )
        .def( init< scitbx::af::const_ref< scitbx::vec3 <double > > const&,
                    scitbx::af::const_ref< int >    const&
                  >
             (( arg_("xyz"),
                arg_("ca_indx")
             ))
            )
        .def("restraint", &w_t::restraint)
        .def("restraint", &w_t::internal_restraint)
        .def("relax", &w_t::relax)
        .def("nmode", &w_t::nmode)
        .def("eigenvalues", &w_t::eigenvalues)
        .def("ca_mode", &w_t::ca_mode)
        .def("mode", &w_t::Project2All)
        .def("project2all", &w_t::project2all)
        .def("Histogram", &w_t::Histogram)
        .def("updateDistArray", &w_t::updateDistArray)
        .def("updateDistArrayAll", &w_t::updateDistArrayAll)
        .def("getDistArray", &w_t::getDistArray)
      ;
    }

  };


  //RTB based ENM //
  struct elastic_rtb_wrapper
  {
    typedef elastic_rtb < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("elastic_rtb", no_init)
        .def( init< scitbx::af::const_ref< scitbx::vec3 <double > > const&,
                    scitbx::af::const_ref< int >    const&,
                    double const&,
                    double const&
                  >
             (( arg_("xyz"),
                arg_("block_start"),
                arg_("cutoff"),
                arg_("scale_factor")
             ))
            )
        .def("nmode", &w_t::nmode)
        .def("mode", &w_t::mode)
      ;
    }

  };
} // namespace <anonymous>

namespace boost_python{

  void wrap_elastic()
  {
   sastbx::refine::elastic_wrapper::wrap();
  }

  void wrap_elastic_rtb()
  {
   sastbx::refine::elastic_rtb_wrapper::wrap();
  }
}

}}
