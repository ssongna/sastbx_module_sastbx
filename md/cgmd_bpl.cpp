#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/module.hpp>

#include <sastbx/md/beads.h>

namespace sastbx { namespace cgmd {

namespace{

  // Residue based CG //
  struct beads_wrapper
  {
    typedef beads < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("beads", no_init)
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
        .def("update", &w_t::update)
        .def("project", &w_t::Project2All)
        .def("restraint", &w_t::restraint)
        .def("relax", &w_t::relax)
      ;
    }

  };


} // namespace <anonymous>

namespace boost_python{

  void wrap_beads()
  {
   sastbx::cgmd::beads_wrapper::wrap();
  }

}

}}
