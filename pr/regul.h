#ifndef SASTBX_REGUL_H
#define SASTBX_REGUL_H

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
#include <boost/math/special_functions/sinc.hpp>
#include <complex>

#include <sastbx/intensity/integrated.h>

using boost::math::sph_bessel;
using scitbx::math::quadrature::gauss_legendre_engine;

namespace sastbx { namespace pr {

  template <typename FloatType>
  class regul_basic
  {
    /*
       * This class implements a regularisation method for p(r) estimation.
       * It allows the input of multiple datasets
       */
    public:
        regul_basic(FloatType const& d_max,
                    int const& n_points,
                    FloatType const& q_max,
                    FloatType const& q_step
                 ):
        d_max_(d_max), n_points_(n_points), width_(d_max/(n_points)), q_max_(q_max), q_step_(q_step), kernels_(d_max, d_max/(n_points), q_max, q_step), n_sets_(0),
        eps_(1e-15)
        {
           SCITBX_ASSERT( n_points_ >4 );
           // make the rbin array please. we'll print bin centers
           for (int ii=0;ii<n_points_;ii++){
             r_.push_back(ii*width_+width_/2.0);
             p_.push_back(0.0);
           }
        }

        /*  Load Data */

        void
        load_data( scitbx::af::const_ref< FloatType > const& this_q,
                   scitbx::af::const_ref< FloatType > const& this_i,
                   scitbx::af::const_ref< FloatType > const& this_s)
        {
          scitbx::af::shared< FloatType > tmp_q;
          scitbx::af::shared< FloatType > tmp_i;
          scitbx::af::shared< FloatType > tmp_s;

          SCITBX_ASSERT( this_q.size() == this_i.size() );
          SCITBX_ASSERT( this_s.size() == this_i.size() );
          for (int ii=0;ii<this_i.size();ii++){
            tmp_q.push_back( this_q[ii] );
            tmp_i.push_back( this_i[ii] );
            SCITBX_ASSERT( this_s[ii]>eps_);
            tmp_s.push_back( this_s[ii]*this_s[ii] ); //make it variances plase
          }
          // and store it locally
          q_.push_back(  tmp_q );
          io_.push_back( tmp_i );
          so_.push_back( tmp_s );

          // also, please update our kernel engine
          kernels_.setup_arrays( this_q );
          scales_.push_back( 1.0 );
          n_sets_++;

        }
        /*  Load the distribution  */
        void
        load_p( scitbx::af::const_ref< FloatType > const& p)
        {
          for (int ii=0;ii<p_.size();ii++){
            p_[ii] = p[ii];
          }
        }

        /*  Load the scales  */
        void
        load_scales( scitbx::af::const_ref< FloatType > const& scales )
        {
          check(scales.size()-1);
          for (int ii=0;ii<scales_.size();ii++){
            scales_[ii]=scales[ii];
          }
        }


        /*  Get r  */
        scitbx::af::shared< FloatType > r()
        {
          return( r_ );
        }

        void check(int const& set_index)
        {
          SCITBX_ASSERT( set_index < n_sets_ );
          SCITBX_ASSERT( set_index >= 0 );
        }

        /* get intensity*/
        scitbx::af::shared< FloatType >
        get_intensity(int const& set_index)
        {
          check(set_index);
          scitbx::af::shared< FloatType > result;
          FloatType scale = scales_[set_index];
          result = kernels_.get_intensity(p_.const_ref(), set_index, scale);
          return( result );
        }

        scitbx::af::shared< FloatType >
        get_delta(int const& set_index)
        {
          check(set_index);
          scitbx::af::shared< FloatType > result;
          scitbx::af::shared< FloatType > delta;
          FloatType scale = scales_[set_index];
          result = kernels_.get_intensity(p_.const_ref(), set_index, scale);
          for (int q_index=0;q_index<result.size();q_index++){
            delta.push_back( (io_[set_index][q_index]-result[q_index])/(so_[set_index][q_index]) );
          }
          return( delta );
        }

        FloatType
        t_int(int const& set_index)
        {
          check(set_index);
          scitbx::af::shared< FloatType > intensity;
          FloatType tmp,result=0;

          intensity = get_intensity(set_index);
          for (int ii=0;ii<intensity.size();ii++){
            tmp = (io_[set_index][ii]-intensity[ii]);
            tmp = tmp*tmp/so_[set_index][ii];
            result += tmp;
          }
          return(result);
        }

        /*  Get derivatives for bins for given index  */
        scitbx::af::shared< FloatType >
        dt_dp( int const& set_index )
        {
           check(set_index);
           scitbx::af::shared< FloatType > result;
           scitbx::af::shared< FloatType > delta;
           FloatType scale = scales_[set_index];
           std::vector< scitbx::af::shared< FloatType > > this_kernel;
           this_kernel = kernels_.get_table_slow_in_r( set_index );
           delta = get_delta(set_index);
           FloatType tmp;

           for (int r_index=0; r_index<this_kernel.size(); r_index++){
             tmp=0;
             for (int q_index=0;q_index<q_[set_index].size();q_index++){
               tmp+=delta[q_index]*this_kernel[r_index][q_index];
             }
             tmp = -tmp*2.0*scale;
             result.push_back(tmp);
           }
           return(result);
        }

        /*  Hessian for Intensity target */
        scitbx::af::shared< FloatType >
        hessian_t( int const& set_index )
        {
          check(set_index);
          FloatType scale = scales_[set_index];
          FloatType tmp;
          int index;
          std::vector< scitbx::af::shared< FloatType >  > this_kernel;
          this_kernel = kernels_.get_table_slow_in_r( set_index );
          scitbx::af::shared< FloatType > result( r_.size()*r_.size(), 0 );
          for (int r_index1=0; r_index1<this_kernel.size(); r_index1++){
            for (int r_index2=r_index1; r_index2<this_kernel.size(); r_index2++){
              tmp=0;
              for (int q_index=0; q_index<q_[set_index].size(); q_index++){
                tmp+=this_kernel[r_index1][q_index] * this_kernel[r_index2][q_index]/so_[set_index][q_index];
              }
              index = r_index1+r_index2*r_.size();
              result[index] = tmp*2.0*scale*scale ;
              index = r_index2+r_index1*r_.size();
              result[index] = tmp*2.0*scale*scale ;
            }
          }
          return( result );
        }

        FloatType
        dt_ds(int const& set_index )
        {
          check(set_index);
          FloatType scale = scales_[set_index];
          FloatType tmp;
          scitbx::af::shared< FloatType > delta;
          scitbx::af::shared< FloatType > intensity;
          delta = get_delta( set_index );
          intensity = kernels_.get_intensity( p_.const_ref(), set_index, 1.0 );
          FloatType result=0;

          for (int ii=0;ii<q_[set_index].size();ii++){
            result += delta[ii]*intensity[ii];
          }
          return ( -2.0*result);
        }

        FloatType
        ddt_dds(int const& set_index )
        {
          check(set_index);
          scitbx::af::shared< FloatType > intensity;
          intensity = kernels_.get_intensity( p_.const_ref(), set_index, 1.0 );
          FloatType result=0;
          for (int ii=0;ii<q_[set_index].size();ii++){
            result += intensity[ii]*intensity[ii]/so_[set_index][ii];
          }
          return ( 2.0*result);
        }


        FloatType
        t_regul()
        {
          FloatType result=0, tmp;
          int n = p_.size()-1;
          result += 0.5*(p_[0]*p_[0]+p_[n]*p_[n]);
          for (int ii=1;ii<p_.size()-1;ii++){
            tmp = p_[ii] - (p_[ii-1]+p_[ii+1])*0.5;
            result += tmp*tmp;
          }
          return(result);
        }


        scitbx::af::shared< FloatType >
        dregul_dp()
        {
          scitbx::af::shared< FloatType > result( p_.size(),0 );

          FloatType tmp_o, tmp_p, tmp_m;
          int n = p_.size()-1;
          // first the end points
          result[0] = p_[0] - (p_[1]   - (p_[0]  +p_[2])*0.5);
          result[n] = p_[n] - (p_[n-1] - (p_[n-2]+p_[n])*0.5);

          // flanking points
          result[1]   = 2.0*(p_[1]   - (p_[0]  +p_[2]  )*0.5) - (p_[2]   - (p_[1]  +p_[3]  )*0.5 ) ;
          result[n-1] = 2.0*(p_[n-1] - (p_[n-2]+p_[n]  )*0.5) - (p_[n-2] - (p_[n-3]+p_[n-1])*0.5 ) ;


          for (int ii=2;ii<n-1;ii++){
            tmp_o = p_[ii]   - (p_[ii-1]+p_[ii+1])*0.5;
            tmp_p = p_[ii-1] - (p_[ii-2]+p_[ii  ])*0.5;
            tmp_m = p_[ii+1] - (p_[ii  ]+p_[ii+2])*0.5;
            result[ii] = 2.0*tmp_o - tmp_m - tmp_p;
          }
          return(result);
        }

        scitbx::af::shared< FloatType >
        hessian_regul()
        {
          if (hessian_regul_.size() == 0){
            scitbx::af::shared< FloatType > result( p_.size()*p_.size(),0 );
            // first the main body please
            for (int ii=2;ii<r_.size()-2;ii++){
              result[ii*r_.size()+ii]      =  3.0;
              result[ii*r_.size()+ii-1]    = -2.0;
              result[ii*r_.size()+ii+1]    = -2.0;
              result[(ii-1)*r_.size()+ii]  = -2.0;
              result[(ii+1)*r_.size()+ii]  = -2.0;
              result[ii*r_.size()+ii-2]    =  0.5;
              result[ii*r_.size()+ii+2]    =  0.5;
              result[(ii-2)*r_.size()+ii]  =  0.5;
              result[(ii+2)*r_.size()+ii]  =  0.5;
            }
            // now the end points
            int n = r_.size();
            result[0]     =  1.5;
            result[1]     = -1.0;
            result[n]     = -1.0;
            result[n+1]   =  2.5;

            result[n*n-1] =  1.5;
            result[n*n-2] = -1.0;

            result[(n-1)*(n) -1] = -1.0;
            result[(n-1)*(n) -2] =  2.5;
            hessian_regul_ = result;
          }
          return(hessian_regul_);
        }



    private:
       // basic variables
       FloatType d_max_, n_points_, q_max_, q_step_, width_;

       // quick integrator
       sastbx::intensity::mset_kernel<FloatType> kernels_;

       // pofr and scales
       scitbx::af::shared< FloatType > p_;
       scitbx::af::shared< FloatType > scales_;

       // user q array
       std::vector< scitbx::af::shared< FloatType > > q_ ;

       // user data
       std::vector< scitbx::af::shared< FloatType > > io_ ;
       std::vector< scitbx::af::shared< FloatType > > so_ ;

       // gradient and snd_der
       std::vector< scitbx::af::shared< FloatType > > d_int_ ;
       std::vector< scitbx::af::shared< FloatType > > H_int_ ;

       scitbx::af::shared< FloatType > d_smooth_;
       scitbx::af::shared< FloatType > H_smooth_;

       scitbx::af::shared< FloatType> r_;
       scitbx::af::shared< FloatType> hessian_regul_;

       int n_sets_;
       FloatType eps_;

  };







}} //namespace sastbx::pr
#endif //SASTBX_regul_H
