#ifndef SASTBX_ZERNIKE_MODEL_H_
#define SASTBX_ZERNIKE_MODEL_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>

#include <scitbx/vec3.h>
#include <scitbx/math/zernike.h>
#include <boost/math/special_functions/bessel.hpp>
#include <complex>

using namespace scitbx::math::zernike;
using boost::math::sph_bessel;

namespace sastbx { namespace zmodel {

  template <typename FloatType>
  class zernike_model
  {
    public:
      zernike_model(
                   scitbx::math::zernike::nlm_array<FloatType> const& C_nlm,
                   scitbx::af::const_ref<FloatType> const& q_array,
                   FloatType const& r_max,
                   int const& n_max
                  ):
                  C_nlm_(C_nlm),
                  C_nn_(n_max),
                  r_max_(r_max),
                  n_max_(n_max),
                  lgf_(n_max)
      {
        nq_=q_array.size();
        for(int i=0;i<nq_;i++)
          q_array_.push_back( q_array[i] );
        build_bnl_array();

      }

     scitbx::af::shared< FloatType >
     calc_intensity_nnl(scitbx::math::zernike::nlm_array<FloatType> nnl_array)
     {
       int start_l;
       scitbx::af::shared<FloatType> I;
       std::complex<FloatType> complex_iq, complex_zero(0,0);
       FloatType tmp_bnl, iq, nnl_coef;
       for(int qq=0; qq<nq_; qq++)
       {
         iq = 0;
         for(int n1=0; n1<=n_max_; n1++)
         {
          for(int n2=0;n2<=n1;n2++)
          { //int nn=(n1<n2)? n1:n2;
           int nn = n2;
           if( is_even( nn ) ) start_l =0;
           else start_l = 1;

           for(int ll=start_l; ll<=nn; ll+=2)
           {
               tmp_bnl= b_n_[n1][qq] * b_n_sign_[n1][ll];
               tmp_bnl=tmp_bnl*b_n_[n2][qq] * b_n_sign_[n2][ll];
               nnl_coef = nnl_array.get_coef(n1,n2,ll).real();
               iq += tmp_bnl*nnl_coef;
           }
          }
         }
        I.push_back(iq);
       }
       return I;
     }


     scitbx::af::shared< FloatType >
     calc_intensity_nlm(scitbx::math::zernike::nlm_array<FloatType> nlm_array)
     {
        calc_invariance_nn( nlm_array );
        return calc_intensity( C_nn_ );
     }

     scitbx::af::shared< FloatType > calc_intensity(scitbx::math::zernike::nl_array<FloatType> nn_array)
     {
       int start_n2;
       scitbx::af::shared<FloatType> I;
       FloatType iq, tmp_bnl;

       for(int qq=0; qq<nq_; qq++)
       {
         iq = 0.0;
         for(int n1=0; n1<=n_max_; n1++)
         {
           if( is_even( n1 ) ) start_n2 =0;
           else start_n2 = 1;

           for(int n2=start_n2; n2<=n1; n2+=2)
           {
             tmp_bnl=b_n_[n1][qq]*b_n_[n2][qq];
             iq += tmp_bnl* nn_array.get_coef(n1,n2);
           }
         }
         I.push_back(iq*2.0); //get the scale right, summed up lower tri-angle and with half values on diag line
       }
       return I;
     }

// utility functions
//
      void calc_invariance_nn(scitbx::math::zernike::nlm_array<FloatType> nlm_array ) {
        FloatType tmp1, tmp2;
        int start_l, start_n2, coef, tmp_n;
        for(int n1=0;n1<=n_max_;n1++) {
          start_n2 = (n1-n1/2*2);
          for(int n2=start_n2; n2<=n1;n2+=2) {
            start_l = (n2-n2/2*2);
            tmp1 = 0;
            for(int l=start_l; l<=n2;l+=2) {
              tmp_n = l-((n1+n2)/2);
              tmp_n = tmp_n-tmp_n/2*2; // tmp_n%2
              coef = pow_1( tmp_n );
              tmp2 = 0;
              for(int m=-l;m<=l;m++) {
                tmp2=tmp2+std::real( std::conj(nlm_array.get_coef(n1,l,m))*nlm_array.get_coef(n2,l,m) );
              }
              tmp1 += tmp2*coef;
            }  //end l
            if( n1 == n2)
              tmp1 = tmp1 / 2.0;  // to scale all pairs avoids double counting in intensity calc
            C_nn_.set_coef(n1,n2,tmp1);
          }  //end n2
        } //end n1
        return;
      }


      void build_bnl_array()
      { // This 3-D array stores coef's b_nl(q) in a way that
        scitbx::af::shared<FloatType> qr_array( nq_, 0.0 );
        FloatType scale;
        for(int qq=0; qq<nq_; qq++)
          qr_array[qq] = q_array_[qq]*r_max_;

        scitbx::af::shared< FloatType > jn( nq_, 0.0 );
        scitbx::af::shared< FloatType > jn1( nq_, 0.0 );
        scitbx::af::shared< FloatType > jn2( nq_, 0.0 );

        for(int qq=0; qq<nq_;qq++) {
          jn[qq]= (sph_bessel(0,qr_array[qq]) );
          jn1[qq]=( sph_bessel(1,qr_array[qq]) );
        }

        for(int nn=0; nn<=n_max_;nn++)
        {
          scale = std::sqrt(2.0*nn+3.0); //scale is due to the zernike radial fn??
          scitbx::af::shared< FloatType > bnq( nq_, 0.0);
          for(int qq=0; qq<nq_; qq++) {
            jn2[qq] = sph_bessel(nn+2,qr_array[qq]);
            bnq[qq] = (jn[qq]+jn2[qq] )/scale;
            jn[qq] = jn1[qq];
            jn1[qq] = jn2[qq];
          }


          b_n_.push_back( bnq );
          scitbx::af::shared< int > b_nl_sign;

          for(int ll=0;ll<=nn;ll+=1)  // to ensure that (nn-ll) is even
          {
            if( is_even( (nn - ll)/2 ) )
              b_nl_sign.push_back(1);
            else
              b_nl_sign.push_back(-1);
          }
          b_n_sign_.push_back( b_nl_sign );
        }
        return;
      }

      inline bool is_even(int num)
      {
        if((num/2)*2 == num) return true;
        else return false;
      }

      inline int pow_1(int num)
      {
        if((num/2)*2 == num) return 1;
        else return -1;
      }


    private:
      scitbx::math::zernike::log_factorial_generator< FloatType > lgf_;
      scitbx::math::zernike::nlm_array<FloatType> C_nlm_;
      scitbx::math::zernike::nl_array<FloatType> C_nn_;
      scitbx::af::shared< scitbx::af::shared< FloatType> > b_n_;
      scitbx::af::shared< scitbx::af::shared<int> > b_n_sign_;
      scitbx::af::shared< FloatType > q_array_;
      FloatType r_max_;
      int l_, n_, n_p_, nq_;
      int n_max_;
  };


}} //namespace sastbx::zmodel
#endif //SASTBX_ZERNIKE_MODEL_H
