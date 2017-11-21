//! Peter Zwart April 05, 2005
#ifndef SASTBX_INTENSITY_DEBYE_H
#define SASTBX_INTENSITY_DEBYE_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <map>

#include <scitbx/vec3.h>
#include <string>
#include <iomanip>

#include <sastbx/intensity/model.h>

namespace sastbx { namespace intensity {

  template <typename FloatType>
  class debye_engine
  {
    public:
      debye_engine(  model <FloatType> const & model,
                     scattering_library< FloatType > const& scat_lib // scattering library
                  ):
      model_(model),
      scat_lib_( scat_lib )
      {
        n_ = model_.size();

        for(int i=1; i < n_; i++) //lower triangle
        { for (int j=0; j < i; j++)
                {
                  distance_matrix_.push_back( model_.distance(i,j) );
                }
        }
      }



// calculate I(s) using debye formula
      scitbx::af::shared<FloatType> I()
      {
        int n_q = scat_lib_.n_q();
        FloatType q, this_I, div,eps=1e-12;
        scitbx::af::shared< FloatType > I_array;
        int accum;

        for (int index=0;index<n_q;index++)
        {
         q  = scat_lib_.q(index);
         sf_i = model_.sf_array_[index][0];
         this_I = sf_i*sf_i/2.0;
         accum=0;
         for(int i=1; i< n_; i++)
         {
          sf_i=model_.sf_array_[index][i];
          this_I += sf_i*sf_i/2.0;
          for(int j=0; j<i; j++)
          {
                sf_j=model_.sf_array_[index][j];
                r = distance_matrix_[accum]; //j+i*(i-1)/2];
                div = r*q;
                if (div>eps){
                this_I += sf_i*sf_j*std::sin(div)/(std::max(div,eps)) ;
                } else {
                this_I += sf_i*sf_j;
                }
                accum++;
          }
         }
         I_array.push_back( 2.0*this_I );
        }

        return (I_array);
      }//end of debye_engine constructor

// retrieve distance between (i,j)
      FloatType Map_Distance(int i, int j)
        {
                return distance_matrix_[j + i*(i-1)/2];
        }

    private:
      model <FloatType> model_;
      scattering_library< FloatType > scat_lib_;
      scitbx::af::shared< FloatType > distance_matrix_;
      int n_, qq;
      FloatType s, I_, sf_i, sf_j, r;
  };





  template <typename FloatType>
  class fast_debye_engine
  {
    public:
      fast_debye_engine(  scitbx::af::const_ref< scitbx::vec3<FloatType> > const& xyz,
                          scitbx::af::const_ref< FloatType > const& q_values,
                          scitbx::af::const_ref< FloatType > const& scattering_factor
                       )
      {
        SCITBX_ASSERT( xyz.size() > 0 );
        SCITBX_ASSERT( q_values.size() > 0 );
        SCITBX_ASSERT( q_values.size() == scattering_factor.size() );
        for (int ii=0;ii<xyz.size();ii++){
          xyz_.push_back( xyz[ii] );
        }
        for (int ii=0;ii<q_values.size();ii++){
          q_.push_back( q_values[ii] );
          sf_sq_.push_back( scattering_factor[ii]*scattering_factor[ii] );
        }
        n_ = xyz.size();
        for(int i=1; i < n_; i++){
          for (int j=0; j < i; j++){
            distance_matrix_.push_back( (xyz_[i]-xyz_[j]).length() );
          }
        }

      }

      FloatType single_i(int this_q_index, FloatType global_b )
      {
        FloatType this_I=0.0, q, div,r;
        FloatType eps=1e-8;
        int accum=0;
        q = q_[this_q_index];

        // don't forget the first atom!
        this_I = 0.5;
        for(int i=1; i< n_; i++){
          this_I += 0.5;
          for(int j=0; j<i; j++){
            r = distance_matrix_[accum];
            div = r*q;
            if (div<eps){
              this_I += 1.0;
            } else {
              this_I += std::sin(div)/(div) ;
            }
            accum++;
          }
        }
        FloatType stol = q_[this_q_index]/(scitbx::constants::pi*4.0);
        this_I=this_I*sf_sq_[this_q_index]*2.0*std::exp(-2.0*global_b*stol*stol);
        return( this_I );
      }


      scitbx::af::shared< FloatType > I(FloatType global_b)
      {
        scitbx::af::shared<FloatType> result;
        for(int ii=0;ii<q_.size();ii++){
          result.push_back( single_i(ii,global_b) );
        }
        return(result);
      }

    private:
      int n_;
      scitbx::af::shared< scitbx::vec3<FloatType> > xyz_;
      scitbx::af::shared< FloatType > q_;
      scitbx::af::shared< FloatType > sf_sq_;
      scitbx::af::shared< FloatType > distance_matrix_;
  };




}}  // namespace sastbx::intensity
#endif // SASTBX_INTENSITY_DEBYE_H
