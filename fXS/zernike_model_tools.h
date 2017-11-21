#ifndef SASTBX_FXS_ZERNIKE_TOOLS_H_
#define SASTBX_FXS_ZERNIKE_TOOLS_H_

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/math/zernike.h>
#include <boost/math/special_functions/bessel.hpp>
#include <complex>
#include <scitbx/wigner/wigner3j.h>


using boost::math::sph_bessel;
using boost::math::legendre_p;
using scitbx::wigner::wigner3j_fast;

namespace array=scitbx::af;

namespace sastbx { namespace fXS {

  template <typename FloatType>
  class zernike_moment_variants
  {
    public:
      zernike_moment_variants() {}
      zernike_moment_variants(
        scitbx::math::zernike::nlm_array<FloatType> const& C_nlm,
        array::const_ref<FloatType> const& q_array,
        FloatType const& r_max,
        int const& n_max,
        int const& l_max
      ): n_max_(n_max), max_L_(l_max), C_nlm_(C_nlm), r_max_(r_max),
         wigner_fast_(n_max*3+1)
      {
        nq_ = q_array.size();
        for(int i=0;i<nq_;i++)
          q_array_.push_back( q_array[i] );
        std::complex<FloatType> complex_I(0,1);
        for(int l=0;l<=n_max_;l++)
          i_pow_.push_back( pow(complex_I,l) );
        setup_b_q();
        compute_G_coef();
      }

// Utility functions
//
// setup_b_q: depending on C_nlm and r_max_
// if any of these two changes, then need to rerun this fn
      void setup_b_q()
      {
        b_n_.clear();
        b_n_sign_.clear();
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
            jn[qq]  = jn1[qq];
            jn1[qq] = jn2[qq];
          }
          b_n_.push_back( bnq );   //to refer: b_n_[n][q_index]

          scitbx::af::shared< int > b_nl_sign;
            //this array goes from 0 to n with increment of 1, to use value, b_n_sign_[n][l]
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

// compute_G_coef: computes Clebsch-Gorden Coefs
// this does not depend on C_nlm nor r_max_
// need to be calculated ONCE only
//
      void compute_G_coef()
      {
        for(int l=0;l<=max_L_;l+=2)
          for(int m=-l;m<=l;m++)
            coef_G_.push_back( compute_G_coef_lm(l,m) );
        return;
      }

      array::shared<FloatType>
      compute_G_coef_lm(int l, int m)  //Clebsch-Gorden Coefs for (l,m)
      {
        FloatType coef_m1, coef_m2, wigner_l1, wigner_l2, wigner_lm1, wigner_lm2, sqrt_value;
        FloatType four_pi=scitbx::constants::pi * 4.0;
        array::shared< FloatType >  coef_lm;  // Clebsch-Gordan Coeffs
        int m2, l2, m1, l1;
        for( l1=0;l1<=n_max_;l1++)
        {
          for( m1 = -l1; m1<=l1; m1++)
          {
            m2 = m1-m; //to satisfy m1-m2-m=0; i.e. m1+m2+m3=0 in standard repr
            for( l2=abs(m2);l2<=n_max_;l2++)
            {
              if( (l1+l2+l) % 2 == 0)
              {
                sqrt_value=sqrt( (2*l1+1)*(2*l2+1)*(2*l+1)/four_pi );
                wigner_l1 = wigner_fast_.compute(l1,l2,l,0,0,0);
                wigner_lm1 = wigner_fast_.compute(l1,l2,l,m1,-m2,-m);
                coef_m1 = wigner_l1* wigner_lm1* sqrt_value;
                coef_lm.push_back( coef_m1 );
                if( m != 0)
                {
                  wigner_l2 = wigner_l1; //wigner_fast_.compute(l2,l1,l,0,0,0);
                  wigner_lm2 = wigner_fast_.compute(l2,l1,l,m2,-m1,-m);
                  coef_m2 = wigner_l2* wigner_lm2* sqrt_value;
                  coef_lm.push_back( coef_m2 );
                }
              }
            }
          }
        }
        return coef_lm;
      }

      array::shared< std::complex<FloatType> >
      compute_a_expansion(int q_index)
      {
        int l,m,n;
        array::shared< std::complex<FloatType> > a_q_lm;
        std::complex< FloatType > tmp_lm_q(0.0,0.0), i_pow_l;
        for(l=0;l<=n_max_;l++)
        {
          i_pow_l = i_pow_[l];
          for(m=-l;m<=l;m++)
          {
            tmp_lm_q = 0.0;
            for(n=l;n<=n_max_;n+=2)  // because n>=l, n,l has the same parity
              tmp_lm_q += b_n_sign_[n][l]*b_n_[n][q_index]*C_nlm_.get_coef(n,l,m)*i_pow_l;
            a_q_lm.push_back( tmp_lm_q );
          }
        }
        return a_q_lm;
      }
//spherical harmonics expansion at q-surface
      array::shared< std::complex<FloatType> >
      compute_i_expansion(int q_index, array::const_ref<std::complex<FloatType> > a_lm)
      {
        int lm_indx(0);
        array::shared< std::complex<FloatType> > I_lm_q;
        for(int l=0;l<=max_L_;l+=2)
          for(int m=-l;m<=l;m++)
          {
            I_lm_q.push_back( compute_I_lm( l, m, a_lm, coef_G_[lm_indx].const_ref() ) );
            lm_indx++;
          }
        return I_lm_q;
      }

// Compute Intensity expansion I_lm //
      std::complex<FloatType>
      compute_I_lm( int l, int m,
      array::const_ref< std::complex<FloatType> > a_q ,
      array::const_ref<FloatType> coef_array)
      {
        //std::cout<<a_q[0]<<" "<<coef_array[0]<<" "<<l<<" "<<m<<std::endl;
        std::complex< FloatType > I_lm_value(0,0), sum_m1, sum_m2;
        FloatType coef_m1, coef_m2, wigner_l1, wigner_l2, wigner_lm1, wigner_lm2;
        int m2,l2, m1, l1, accum(0);


        for( l1=0;l1<=n_max_;l1++)
        {
          for( m1 = -l1; m1<=l1; m1++)
          {
            m2 = m1-m; //to satisfy m1-m2-m=0; i.e. m1+m2+m3=0 in standard repr
            for( l2=abs(m2);l2<=n_max_;l2++)
            {
              if( is_even(l1+l2+l) )
              {
                coef_m1 = coef_array[accum++];
//                if(coef_m1==0) continue;
                sum_m1 = coef_m1*a_q[l1*l1+l1+m1]*conj( a_q[l2*l2+l2+m2] );
                if( is_even(m1) )
                  I_lm_value += sum_m1;
                else
                  I_lm_value -= sum_m1;

                if( m != 0)
                {
                  coef_m2 = coef_array[accum++];
                  sum_m2 = coef_m2*a_q[l2*l2+l2+m2]*conj( a_q[l1*l1+l1+m1] );
                  if( is_even(m2) )
                    I_lm_value += sum_m2;
                  else
                    I_lm_value -= sum_m2;
                }
              }
            }
          }
        }

        if( is_even(m) )
          return I_lm_value;
        else
          return -I_lm_value;
      }

      array::shared< std::complex<FloatType > > get_a_expansion( int q_index )
      {
        return compute_a_expansion(q_index);
      }

      array::shared< std::complex<FloatType > >
      get_i_expansion( int q_index, array::const_ref<std::complex<FloatType> > a_lm )
      {
        return compute_i_expansion(q_index, a_lm);
      }

      array::shared< FloatType >
      get_all_blq(scitbx::math::zernike::nlm_array<FloatType> & C_nlm)
      {
         C_nlm_.load_coefs( C_nlm.nlm(), C_nlm.coefs().const_ref() );
         return get_all_blq2();
      }
      array::shared< FloatType > get_all_blq2()
      {
         array::shared< FloatType > B;
         array::shared< FloatType > blq;
         for(int i=0;i<q_array_.size();i++)
         {
           blq = get_real_coef(i);
           B.insert( B.end(), blq.begin(), blq.end() );
         }
         return B;
      }

      array::shared< FloatType > get_real_coef(int q_index) //B_ls
      {
        array::shared< FloatType > B_l_q;
        array::shared< std::complex< FloatType > > a_lm;
        array::shared< std::complex< FloatType > > i_lm;
        FloatType this_b;
        int indx(0);

        a_lm = compute_a_expansion(q_index);
        i_lm = compute_i_expansion(q_index, a_lm.const_ref() );
        for(int l=0;l<=max_L_;l+=2)
        {
          this_b = 0.0;
          for(int m=-l;m<=l;m++)
            this_b += norm(i_lm[indx++]);
          B_l_q.push_back(this_b);
        }
        //std::cout<<B_l_q[0]<<" "<<a_lm[0]<<" "<<i_lm[0]<<" "<<q_array_[q_index]<<std::endl;
        return B_l_q;
      }

      inline bool is_even(int n)  { return (n==n/2*2); }

    private:
      scitbx::math::zernike::nlm_array<FloatType> C_nlm_;
      scitbx::af::shared< scitbx::af::shared< std::complex<FloatType> > > a_; //a_lm array
      scitbx::af::shared< scitbx::af::shared< std::complex<FloatType> > > i_; //i_lm array
      scitbx::af::shared< array::shared< FloatType > > coef_G_;  // Clebsch-Gordan Coeffs
      scitbx::af::shared< scitbx::af::shared< FloatType > > B_; // expansion coefs: B_l array
      scitbx::af::shared< scitbx::af::shared< FloatType> > b_n_;  //zernike radial integration
      scitbx::af::shared< scitbx::af::shared<int> > b_n_sign_;  // sign given by (-1)^((n-l)/2)
      array::shared<FloatType> q_array_;
      array::shared<std::complex<FloatType> > i_pow_;
      wigner3j_fast<FloatType> wigner_fast_;
      FloatType r_max_;
      int n_max_, nq_, max_L_;
  };

}} //namespace sastbx::fXS
#endif //SASTBX_FXS_ZERNIKE_TOOLS_H_
