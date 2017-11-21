#ifndef SASTBX_RBENGINE_H
#define SASTBX_RBENGINE_H
#endif

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <map>

#include <string>
#include <iomanip>
#include <scitbx/vec3.h>
#include <sastbx/rigid_body/rigid_body.h>


namespace sastbx { namespace rigid_body {

   template<typename FloatType>
   class rb_engine //It should be capable to handle k subunits
   {
        public:
	  rb_engine( scitbx::af::const_ref< rigid_body< FloatType> > const& rbs,
		     int const& nbody): nbody_(nbody)
          {
	    SCITBX_ASSERT(rbs.size() == nbody_);

	  }


	  scitbx::af::shared <FloatType> get_pr()
	  {
	  }

 
	  void get_intensity()
	  {
	  }

          
        private:
	  int nbody_;
	  scitbx::af::shared< FloatType > histogram;

   }; //end of intensity class



} //namepsace rigid_body
} //namespace sastbx

