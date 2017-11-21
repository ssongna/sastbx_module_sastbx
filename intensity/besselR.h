// ! Haiguang Liu Mar 11, 2009
#ifndef SASTBX_BESSEL_R_H
#define SASTBX_BESSEL_R_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <map>

#include <scitbx/vec3.h>
#include <string>
#include <iomanip>
#include <scitbx/math/quadrature.h>
#include <boost/math/special_functions/bessel.hpp>
#include <complex>

using boost::math::sph_bessel;

namespace sastbx { namespace intensity {

  template <typename FloatType>
  class besselr
  {
    public:
        besselr(FloatType const& r_step,
                   FloatType const& r_max,
                   int const& l
                 ):
        r_step_(r_step),r_max_(r_max),l_(l)
        {
        FloatType x,half_step,f_i;
        half_step = r_step_/2.0;

        for(int i=0;i<=l_;i++)
        {
          scitbx::af::shared <FloatType> bessel_l;
          for(x=0;x<1.0;x+=half_step)
          { bessel_l.push_back( sph_bessel(i,x) ); }
          for(x=x;x<=r_max_;x+=r_step)
          {
                bessel_l.push_back( sph_bessel(i,x) );
          }
          value_.push_back(bessel_l);
        }
        }

    FloatType get_value(int i, FloatType r)
    {
        int indx;
        FloatType result,r0,r1,step;
        if(r<1.0)
        {step=r_step_/2.0;
        indx = static_cast<int> (r/step);
        r1 = static_cast<FloatType>(step*indx);
        result = value_[i][indx]+(value_[i][indx+1]-value_[i][indx])*(r1+step-r)/step;
        return result;
        }
        else {step=r_step_;
        indx = static_cast<int> ((1.0+r)/step);
        r1 = static_cast<FloatType>(step*indx)-1.0;
        result = value_[i][indx]+(value_[i][indx+1]-value_[i][indx])*(r1-r+step)/step;
        return result;
        }
    }
    void print_out()
    {
        FloatType x;
        FloatType step=r_step_;
        int indx;
        for(int i=0;i<=l_;i++)
        {
         indx=static_cast<int>(1.0/step*2);
         for(x=1.0;x<=r_max_;x+=step)
         {
//              std::cout<<x<<" "<<get_value(i,x-0.05)<<" "<<sph_bessel(i,x-0.05)<<std::endl;
                std::cout<<x<<"\t"<<value_[i][indx]<<" "<<sph_bessel(i,x)<<std::endl;
                indx++;
         }
        }
        return;
    }
    private:
        FloatType r_max_, r_step_;
        scitbx::af::shared < scitbx::af::shared <FloatType> > value_;
        int l_;
  };
}} //namespace sastbx::intensity
#endif //SASTBX_BESSEL_R_H
