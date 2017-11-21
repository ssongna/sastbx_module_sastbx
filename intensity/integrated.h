// ! Haiguang Liu Mar 11, 2009
#ifndef SASTBX_INTEGRATED_H
#define SASTBX_INTEGRATED_H

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

using boost::math::sph_bessel;
using scitbx::math::quadrature::gauss_legendre_engine;

namespace sastbx { namespace intensity {

  template <typename FloatType>
  class integrated
  {
    public:
        integrated(FloatType const& f_min,
                   FloatType const& f_max,
                   FloatType const& delta,
                   FloatType const& f_step,
                   FloatType const& q_step,
                   FloatType const& q_max,
                   int const& l
                 ):
        f_min_(f_min), f_max_(f_max), q_max_(q_max), q_step_(q_step), delta_(delta),resl_(f_step), l_(l), n_p_(10),np_f_((f_max-f_min)/f_step),np_q_(q_max/q_step+1.0)
        {
/*
 * for f in [f_min, f_max]
 * integrate J_l(r*q)*r^2, from f to f+delta
 * q is in a predefined range and resolution
 */
        gauss_legendre_engine <FloatType> gle(n_p_);
        w_=gle.w();
        x_=gle.x();

        FloatType tmp_q,  tmp_f, sum, a, b, coef,off, r;
        int accum;
        for(tmp_f=f_min_;tmp_f<=f_max_+resl_; tmp_f += resl_ )
        {
        scitbx::af::shared <FloatType> Inte_q;
         a= tmp_f;
         b= tmp_f + delta_;
         coef=(b-a) / 2.0;
         off =(b+a) / 2.0;
         for(tmp_q=0.00;tmp_q<=(q_max_+q_step_*2); tmp_q+=q_step_)
         {
          sum = 0.0;
          for(int i=0;i<n_p_;i++)
          {
                r=coef*x_[i] + off;
                sum += w_[i]*(sph_bessel(l,tmp_q*r)*r*r);
          }
          sum *=coef;
          Inte_q.push_back(sum);
         }//end tmp_q
        Integral_.push_back(Inte_q);
        } //end tmp_f
        }

        FloatType Get_integral(FloatType const& f_w, FloatType const& q)
        {
         FloatType result,f1,f2,q1,q2;
         SCITBX_ASSERT( q>=0.0 && q <= q_max_);
         SCITBX_ASSERT( f_w <= f_max_ && f_w >= f_min_);
         int indx_q = static_cast<int>(q/q_step_);
         int indx_f = static_cast<int>((f_w-f_min_)/resl_);
         if(indx_q == (q/q_step_))
         {
                if(indx_f == (f_w-f_min_) /resl_)
                { return Integral_[indx_f][indx_q];}
                else
                {
                  f1=indx_f*resl_+f_min_; f2=f1+resl_;
                  return (Integral_[indx_f][indx_q]*(f_w-f1) + Integral_[indx_f+1][indx_q]*(f2-f_w))/resl_;
                }
         }
         result = interpolate(indx_f, indx_q, f_w-f_min_, q);
         return result;
        }

        FloatType interpolate(int const& indx_x, int const& indx_y, FloatType const& x, FloatType const& y)
        {
         FloatType x1,x2,y1,y2, coef, result;
         FloatType f11,f12,f21,f22;
         x1 = static_cast<FloatType> (indx_x) *resl_;
         x2 = x1 + resl_;
         y1 = static_cast<FloatType> (indx_y) *q_step_;
         y2 = y1 + q_step_;
         f11 = Integral_[indx_x][indx_y];
         f12 = Integral_[indx_x][indx_y+1];
         f21 = Integral_[indx_x+1][indx_y];
         f22 = Integral_[indx_x+1][indx_y+1];
         coef = 1.0/q_step_/resl_;
         result = f11*(x2-x)*(y2-y) + f21*(x-x1)*(y2-y) + f12*(x2-x)*(y-y1) + f22*(x-x1)*(y-y1);
         result *=coef;

         return result;
        }
        void print_out()
        {
         FloatType tmp_f, tmp_q;
         int index_f, index_q;
         for(tmp_f=f_min_;tmp_f<f_max_; tmp_f += resl_ )
         {
                index_f=static_cast<int> ((tmp_f-f_min_)/resl_);
                for(index_q=0;index_q<np_q_;index_q++)
                {
                std::cout<<std::setw(12)<<Integral_[index_f][index_q]<<"\t";
                std::cout<<std::setw(12)<<Get_integral(tmp_f,index_q*q_step_)<<"\n";
                }
         }
        return ;
        }
    private:
        scitbx::af::shared< scitbx::af::shared<FloatType> > Integral_;
        FloatType f_min_, f_max_, q_max_,resl_,delta_, q_step_;
        int l_, n_p_;
        FloatType np_f_,np_q_;
        scitbx::af::shared< FloatType > w_, x_;
  };


  template <typename FloatType>
  class integrated_sinc
  {
    public:
        integrated_sinc(FloatType const& f_min,
                        FloatType const& f_max,
                        FloatType const& delta,
                        FloatType const& f_step,
                        FloatType const& q_step,
                        FloatType const& q_max
                 ):
        f_min_(f_min), f_max_(f_max), q_max_(q_max), q_step_(q_step), delta_(delta), resl_(f_step), n_p_(40),np_f_((f_max-f_min)/f_step),np_q_(q_max/q_step+1.0)
        {
        gauss_legendre_engine <FloatType> gle(n_p_);
        w_=gle.w();
        x_=gle.x();

        FloatType tmp_q,  tmp_f, sum, a, b, coef,off, r;
        int accum;
        for(tmp_f=f_min_;tmp_f<=f_max_+resl_; tmp_f += resl_ )
        {
        scitbx::af::shared <FloatType> Inte_q;
     //   Inte_q.reserve(static_cast<int>(np_q_+3));
         a= tmp_f;
         b= tmp_f + delta_;
         coef=(b-a) / 2.0;
         off =(b+a) / 2.0;
         for(tmp_q=0.00;tmp_q<=(q_max_+q_step_*2); tmp_q+=q_step_)
         {
          sum = 0.0;
          for(int i=0;i<n_p_;i++)
          {
                r=coef*x_[i] + off;
                sum += w_[i]*( boost::math::sinc_pi<FloatType>(tmp_q*r) );
          }
          sum *=coef;
	  //TESTING for the scaling, see Forster 2008, JMB, summplementary
	  sum *= exp(-0.23*tmp_q*tmp_q);
          Inte_q.push_back(sum);
         }//end tmp_q
        Integral_.push_back(Inte_q);
        } //end tmp_f
        }


        FloatType Get_integral(FloatType const& f_w, FloatType const& q)
        {
         FloatType result,f1,f2,q1,q2;
         SCITBX_ASSERT( q>=0.0 && q <= q_max_);
         SCITBX_ASSERT( f_w <= f_max_ && f_w >= f_min_);
         int indx_q = static_cast<int>(q/q_step_);
         int indx_f = static_cast<int>((f_w-f_min_)/resl_);
         if(indx_q == (q/q_step_))
         {
                if(indx_f == (f_w-f_min_) /resl_)
                { return Integral_[indx_f][indx_q];}
                else
                {
                  f1=indx_f*resl_+f_min_; f2=f1+resl_;
                  return (Integral_[indx_f][indx_q]*(f_w-f1) + Integral_[indx_f+1][indx_q]*(f2-f_w))/resl_;
                }
         }
         result = interpolate(indx_f, indx_q, f_w-f_min_, q);
         return result;
        }

        FloatType interpolate(int const& indx_x, int const& indx_y, FloatType const& x, FloatType const& y)
        {
         FloatType x1,x2,y1,y2, coef, result;
         FloatType f11,f12,f21,f22;
         x1 = static_cast<FloatType> (indx_x) *resl_;
         x2 = x1 + resl_;
         y1 = static_cast<FloatType> (indx_y) *q_step_;
         y2 = y1 + q_step_;
         f11 = Integral_[indx_x][indx_y];
         f12 = Integral_[indx_x][indx_y+1];
         f21 = Integral_[indx_x+1][indx_y];
         f22 = Integral_[indx_x+1][indx_y+1];
         coef = 1.0/q_step_/resl_;
         result = f11*(x2-x)*(y2-y) + f21*(x-x1)*(y2-y) + f12*(x2-x)*(y-y1) + f22*(x-x1)*(y-y1);
         result *=coef;

         return result;
        }
        void print_out()
        {
         FloatType tmp_f, tmp_q;
         int index_f, index_q;
         for(tmp_f=f_min_;tmp_f<f_max_; tmp_f += resl_ )
         {
                index_f=static_cast<int> ((tmp_f-f_min_)/resl_);
                for(index_q=0;index_q<np_q_;index_q++)
                {
                std::cout<<std::setw(12)<<Integral_[index_f][index_q]<<"\t";
                std::cout<<std::setw(12)<<Get_integral(tmp_f,index_q*q_step_)<<"\n";
                }
         }
        return ;
        }
    private:
        scitbx::af::shared< scitbx::af::shared<FloatType> > Integral_;
        FloatType f_min_, f_max_, q_max_,resl_,delta_, q_step_;
        int l_, n_p_;
        FloatType np_f_,np_q_;
        scitbx::af::shared< FloatType > w_, x_;
  };


  template <typename FloatType>
  class block_integrator
  {
    public:
        block_integrator( FloatType const& dmax,
                          FloatType const& width,
                          FloatType const& q_max,
                          FloatType const& q_step ):
        dmax_(dmax),
        width_(width),
        q_max_(q_max),
        q_step_(q_step),
        integrator_(0,dmax,width,width/5.0,q_step,q_max) // the term width/5.0 is needed to provide sufficient numeric accuracy.
                                                         // These settings are good for a sphere with diameter 250 and q values up to 1 (!)
        {
         SCITBX_ASSERT(q_step>0);
         // setup the q array used for interpolation the final result from
         int q_size = int(q_max/q_step) + 1;
         for (int ii=0;ii<q_size;ii++){
           q_.push_back( ii*q_step );
         }

         // make the rbin array please
         int n=int(dmax_/width_+0.5);
         for (int ii=0;ii<n;ii++){
          r_.push_back(ii*width_);
         }

         // make the partial intensities please
         for (int ii=0;ii<q_.size();ii++){
           scitbx::af::shared<FloatType> tmp;
           for (int jj=0;jj<n;jj++){
              tmp.push_back(integrator_.Get_integral(r_[jj], q_[ii]));
           }
           partial_i_.push_back( tmp );
         }

         // lets make a transposed array as well. use for gradients later
         for (int jj=0;jj<n;jj++){
           scitbx::af::shared<FloatType> tmp;
           partial_i_flipped_.push_back( tmp );
           for (int ii=0;ii<q_.size();ii++){
             partial_i_flipped_[jj].push_back( partial_i_[ii][jj] );
           }
         }
         // all done
        }

    // here we setup the arrays
    void
    setup_arrays( scitbx::af::const_ref<FloatType> const& user_q_array )
    {

       // here we setup the interpolation aid
       scitbx::af::shared< FloatType > tmp_q;
       scitbx::af::shared< int > tmp_bin;
       scitbx::af::shared< FloatType > tmp_p;
       int n;
       for (int ii=0;ii<user_q_array.size();ii++){
         tmp_q.push_back( user_q_array[ii] );
         n = int(user_q_array[ii]/q_step_); // we round off to the lowest nearest integer (floor)
         tmp_bin.push_back( n );
         tmp_p.push_back( (user_q_array[ii]-q_[n])/q_step_ );
       }
       user_q_array_ = tmp_q;
       user_q_bin_ = tmp_bin;
       p_array_ = tmp_p;
    }


    FloatType
    three_point_interpolate( int const& q_index )
    {
      FloatType fmo,f,fpo,p, result;
      p = p_array_[q_index] ; // user_q_array_[q_index]-q_[user_q_bin_[q_index]];
      f = intensity_array_[user_q_bin_[q_index]];
      fmo = intensity_array_[user_q_bin_[q_index]-1];
      fpo = intensity_array_[user_q_bin_[q_index]+1];
      result = ((p-1.0)*fmo+(p+1.0)*fpo)*0.5*p  + (1.0-p*p)*f;
      return(result);
    }

    FloatType
    two_point_interpolate( int const& q_index ) // forward interpolation always works because we do a floor operator
    {
      FloatType f,fpo,p, result;
      p = user_q_array_[q_index]-q_[user_q_bin_[q_index]];
      p = p/q_step_;
      f = intensity_array_[user_q_bin_[q_index]];
      fpo = intensity_array_[user_q_bin_[q_index]+1];
      result = (1.0-p)*f + p*fpo;
      return(result);
    }

    FloatType
    interpolate(int q_index)
    {
      FloatType result=0;
      if (q_index > 0){
        result = three_point_interpolate( q_index );
      } else {
        result = two_point_interpolate( q_index );
      }
      return( result );
    }

    scitbx::af::shared< FloatType>
    get_intensity(scitbx::af::const_ref<FloatType> const& weights)
    {
      //---------------------------
      //make sure all is arrays size wise
      SCITBX_ASSERT( weights.size() == r_.size() );

      //---------------------------
      // we first make a lookup table for
      // intensities used in interpolation
      scitbx::af::shared<FloatType> result_for_interpolation;
      FloatType tmp_result;
      for (int ii=0;ii<q_.size();ii++){
        tmp_result=0;
        for (int jj=0;jj<r_.size();jj++){
          tmp_result+= weights[jj]*partial_i_[ii][jj]; // sum the partial intensities
        }
        result_for_interpolation.push_back( tmp_result );
      }
      intensity_array_ = result_for_interpolation;
      //---------------------------
      // Now we need to get intensities via interpolation
      scitbx::af::shared< FloatType > result;
      for (int ii=0;ii<user_q_array_.size();ii++){
        tmp_result = interpolate( ii );
        result.push_back( tmp_result );
      }
      return( result );
    }

    private:
       integrated_sinc<FloatType> integrator_;
       FloatType dmax_, width_, q_step_, q_max_;
       scitbx::af::shared< FloatType > q_;
       scitbx::af::shared< FloatType > r_;
       scitbx::af::shared< FloatType > intensity_array_;
       std::vector< scitbx::af::shared< FloatType> > partial_i_; // partial_i_[ q_indices ][ weight_indices ]
       std::vector< scitbx::af::shared< FloatType> > partial_i_flipped_; // partial_i_flipped_[ weight_indices ][ q_indices ]

       //
       scitbx::af::shared< int > user_q_bin_;
       scitbx::af::shared< FloatType> user_q_array_;
       scitbx::af::shared< FloatType> p_array_;


  };


  // This class is as the above one, but it allows one to
  // load multiple q ranges, possibly with different q values
  // This is to provide support for multi set integrator with one
  // specific set of kernels
  template <typename FloatType>
  class mset_kernel
  {
    public:
        mset_kernel( FloatType const& dmax,
                               FloatType const& width,
                               FloatType const& q_max,
                               FloatType const& q_step ):
        dmax_(dmax),
        width_(width),
        q_max_(q_max),
        q_step_(q_step),
        n_sets_loaded_(0),
        integrator_(0,dmax+width,width,width/20.0,q_step,q_max)
        {
         SCITBX_ASSERT(q_step>0);
         // setup the q array used for interpolation the final result from
         int q_size = int(q_max/q_step) + 1;
         for (int ii=0;ii<q_size;ii++){
           q_.push_back( ii*q_step );
         }
         // make the rbin array please
         int n=int(dmax_/width_+0.5);
         for (int ii=0;ii<n;ii++){
          r_.push_back(ii*width_);
         }

         // make the partial intensities please
         for (int ii=0;ii<q_.size();ii++){
           scitbx::af::shared<FloatType> tmp;
           for (int jj=0;jj<n;jj++){
              tmp.push_back(integrator_.Get_integral(r_[jj], q_[ii]));
           }
           partial_i_.push_back( tmp );
         }

         // lets make a transposed array as well. use for gradients later
         for (int jj=0;jj<n;jj++){
           scitbx::af::shared<FloatType> tmp;
           partial_i_flipped_.push_back( tmp );
           for (int ii=0;ii<q_.size();ii++){
             partial_i_flipped_[jj].push_back( partial_i_[ii][jj] );
           }
         }
         // all done
        }



    // here we setup the arrays
    void
    setup_arrays( scitbx::af::const_ref<FloatType> const& user_q_array )
    {
       n_sets_loaded_++;
       // here we setup the interpolation aid
       scitbx::af::shared< FloatType > tmp_q;
       scitbx::af::shared< int > tmp_bin;
       scitbx::af::shared< FloatType > tmp_p;
       int n;
       for (int ii=0;ii<user_q_array.size();ii++){
         tmp_q.push_back( user_q_array[ii] );
         n = int(user_q_array[ii]/q_step_); // we round off to the lowest nearest integer (floor)
         tmp_bin.push_back( n );
         tmp_p.push_back( (user_q_array[ii]-q_[n])/q_step_ );
       }
       user_q_array_.push_back( tmp_q );
       user_q_bin_.push_back( tmp_bin );
       p_array_.push_back( tmp_p );

       // now compute the tables
       precompute_partial_intensities( n_sets_loaded_-1 );

    }




    //--------------------------------------------
    FloatType
    three_point_interpolate( int const& q_index, int const& set_index )
    {
      FloatType fmo,f,fpo,p, result;
      p = p_array_[set_index][q_index];
      f = intensity_array_[user_q_bin_[set_index][q_index]];
      fmo = intensity_array_[user_q_bin_[set_index][q_index]-1];
      fpo = intensity_array_[user_q_bin_[set_index][q_index]+1];
      result = ((p-1.0)*fmo+(p+1.0)*fpo)*0.5*p  + (1.0-p*p)*f;
      return(result);
    }

    FloatType
    two_point_interpolate( int const& q_index, int const& set_index) // forward interpolation always works because we do a floor operator
    {
      FloatType f,fpo,p, result;
      p = user_q_array_[set_index][q_index]-q_[user_q_bin_[set_index][q_index]];
      p = p/q_step_;
      f = intensity_array_[ user_q_bin_[set_index][q_index] ];
      fpo = intensity_array_[ user_q_bin_[set_index][q_index]+1 ];
      result = (1.0-p)*f + p*fpo;
      return(result);
    }

    FloatType
    interpolate(int const& q_index, int const& set_index)
    {
      FloatType result=0;
      if (q_index > 2){
        result = three_point_interpolate( q_index, set_index );
      } else {
        result = two_point_interpolate( q_index, set_index);
      }
      return( result );
    }
    //--------------------------------------------




    void precompute_partial_intensities( int const& set_index )
    {
       // here we interpolate out the table we have to the user supplied q values
       FloatType tmp_result;
       std::vector< scitbx::af::shared<FloatType> > this_table;

       for (int jj=0;jj<r_.size();jj++){

         // first we assign the proper partial intensities to the 'shared' intensity array used in the interpolation routines
         intensity_array_ = partial_i_flipped_[ jj ]; //

         // now please carry out the interpolation for this value of r
         scitbx::af::shared<FloatType> partial_i_for_given_r_and_set_on_user_q_range;
         for (int ii=0;ii<user_q_array_[set_index].size();ii++){
           tmp_result = interpolate( ii, set_index );
           partial_i_for_given_r_and_set_on_user_q_range.push_back( tmp_result );
         }

         // push this guy onto the table
         this_table.push_back(partial_i_for_given_r_and_set_on_user_q_range);

       }

       // Now the table has been made for this set, store it please
       partial_i_flipped_interpolated_.push_back( this_table );
       unflip( set_index );

    }

    void unflip(int const& set_index)
    {
      std::vector< scitbx::af::shared<FloatType> > this_table;
      for (int ii=0;ii<user_q_array_[set_index].size();ii++){
        scitbx::af::shared<FloatType> tmp( r_.size(), 0 );
        this_table.push_back( tmp );
      }
      // table has been build, now fill up the entries please
      for (int jj=0;jj<r_.size();jj++){
        for (int ii=0;ii<user_q_array_[set_index].size();ii++){
           this_table[ii][jj]= partial_i_flipped_interpolated_[set_index][jj][ii];
        }
      }
      partial_i_interpolated_.push_back( this_table );
    }

    // This function allows us to compute an intensity curve of a given p(r)
    scitbx::af::shared< FloatType>
    get_intensity(scitbx::af::const_ref<FloatType> const& weights, int const& set_index, FloatType const& scale)
    {
      //---------------------------
      //make sure all is good arrays size wise
      SCITBX_ASSERT( weights.size() == r_.size() );
      // we now have to sum the partial intensities per q values
      scitbx::af::shared<FloatType> result;
      for (int ii=0;ii<user_q_array_[set_index].size();ii++){
        FloatType tmp=0;
        // we have to loop over r
        for (int jj=0;jj<r_.size();jj++){
          tmp += partial_i_interpolated_[set_index][ii][jj]*weights[jj]; // slow in q
        }
        result.push_back( tmp*scale );
        //result.push_back( tmp*scale*exp(-0.23*user_q_array_[set_index][ii]^2.0) );
      }
      return( result );
    }


    // this function returns the appropriate lookup table
    std::vector<  scitbx::af::shared< FloatType > >
    get_table_slow_in_q(int const& set_index )
    {
      return(partial_i_interpolated_[set_index]);

    }

    // this function returns the appropriate lookup table
    std::vector<  scitbx::af::shared< FloatType > >
    get_table_slow_in_r(int const& set_index )
    {
      return(partial_i_flipped_interpolated_[set_index]);

    }



    private:
       integrated_sinc<FloatType> integrator_;
       FloatType dmax_, width_, q_step_, q_max_;
       scitbx::af::shared< FloatType > q_;
       scitbx::af::shared< FloatType > r_;
       scitbx::af::shared< FloatType > intensity_array_;
       std::vector< scitbx::af::shared< FloatType> > partial_i_;         // partial_i_[ q_indices ][ weight_indices ]
       std::vector< scitbx::af::shared< FloatType> > partial_i_flipped_; // partial_i_flipped_[ weight_indices ][ q_indices ]


       std::vector< std::vector< scitbx::af::shared< FloatType> > > partial_i_interpolated_; // precomputed intepolated kernel values; slow in q
       std::vector< std::vector< scitbx::af::shared< FloatType> > > partial_i_flipped_interpolated_;// precomputed intepolated kernel values, flipped; slow in r

       int n_sets_loaded_;
       std::vector< scitbx::af::shared< int > > user_q_bin_;
       std::vector< scitbx::af::shared< FloatType> > user_q_array_;
       std::vector< scitbx::af::shared< FloatType> > p_array_;


  };








}} //namespace sastbx::intensity
#endif //SASTBX_INTEGRATED_H
