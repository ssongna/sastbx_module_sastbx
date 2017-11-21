//! Peter Zwart April 05, 2005
#ifndef SASTBX_DATA_REDUCTION_H
#define SASTBX_DATA_REDUCTION_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <cctbx/xray/scatterer_flags.h>
#include <cctbx/sgtbx/site_symmetry.h>
#include <cctbx/adptbx.h>
#include <cctbx/xray/scatterer.h>

#include <map>

#include <scitbx/vec2.h>
#include <iostream>
#include <iomanip>

/*
 *  Basic assumption: x is the slow index, y the fast
 */


namespace sastbx { namespace data_reduction {
  template <typename FloatType>
  class perpendicular_detector
  {
    public:
      perpendicular_detector( int const& nx, int const& ny,
                              FloatType const& scalex, FloatType const& scaley,
                              FloatType const& orgx, FloatType const& orgy,
                              FloatType const& distance, FloatType const& wavelength
                              )
      {
         nx_ = nx;
         ny_ = ny;
         scalex_ = scalex;
         scaley_ = scaley;
         orgx_ = orgx;
         orgy_ = orgy;
         distance_ = distance; // in mm please
         wavelength_ = wavelength; // in angstrom please

         // set binning related variables to -9 (for pickling getstate)
         start_q_ = -9;
         end_q_   = -9;
         q_step_  = -9;

         SCITBX_ASSERT( nx > 0 );
         SCITBX_ASSERT( ny > 0 );
         min_q_ = 1e5;
         max_q_=-1;
         FloatType tmp1,tmp2,result;
         for (int ii=0;ii<nx_;ii++){
           tmp1 = (ii-orgx_)*scalex_;
           for (int jj=0;jj<ny_;jj++){
             tmp2 = (jj-orgy_)*scaley_;
             // first we need to compute the angle
             result = std::sqrt( tmp1*tmp1 + tmp2*tmp2 );
             result = std::atan2( result, distance_ )/2.0;
             // now compute q
             result = scitbx::constants::pi*4.0*std::sin(result)/wavelength_;
             if (result > max_q_){
               max_q_= result;
             }
             if (result < min_q_){
               min_q_ = result;
             }
             q_.push_back( result );
           }
         }

      }

    int nx(){ return(nx_); }
    int ny(){ return(ny_); }
    FloatType scale_x(){ return(scalex_); }
    FloatType scale_y(){ return(scaley_); }
    FloatType origin_x(){ return(orgx_); }
    FloatType origin_y(){ return(orgy_); }
    FloatType distance(){ return(distance_); }
    FloatType wavelength(){ return(wavelength_); }
    scitbx::af::shared<FloatType> q(){ return(q_); }

    void compute_bin_indices(FloatType q_step, FloatType start_q, FloatType end_q)
    {
      // declare some bin related data structures
      n_bins_ = static_cast<int>( std::floor(1+(end_q - start_q)/q_step) );
      q_step_ = q_step;
      start_q_ = start_q;
      end_q_ = end_q;
      for (int ii=0;ii<=n_bins_;ii++){
        std::vector<int> tmp;
        q_low_bin_indices_.push_back( tmp );
        q_low_bin_values_.push_back( start_q_ + q_step_*ii );
        q_mean_bin_values_.push_back(0.0);
      }
      // now do the actual binning
      FloatType this_q, ft_q;
      int i_q;
      for (int ii=0;ii<q_.size();ii++){
        this_q = q_[ ii ];
        // figure out in which bin this guy falls
        ft_q = (this_q-start_q_)/q_step_;
        i_q  = std::floor( ft_q );
        if (i_q <=0 ){ // if q falls before our limit
          before_indices_.push_back( ii );
        }
        if (i_q>n_bins_){ // if it falls after our specificied range
          after_indices_.push_back( ii );
        }
        if (i_q >= 0) {  if (i_q<=n_bins_) {
          q_low_bin_indices_[ i_q ].push_back( ii ); // push back the bin number please
          q_mean_bin_values_[ i_q ] += this_q;
        } }
      }

      // get means
      for (int ii=0;ii<q_mean_bin_values_.size();ii++){
        q_mean_bin_values_[ii]/=(q_low_bin_indices_[ii].size()+1e-13);;
      }
      // all done
    }

    scitbx::af::shared< FloatType > q_low_bin_values()
    {
      return( q_low_bin_values_ );
    }

    scitbx::af::shared< FloatType > q_mean_bin_values()
    {
      return( q_mean_bin_values_);
    }

    scitbx::af::shared< std::vector<int> > q_bin_indices()
    {
      return( q_low_bin_indices_ );
    }

    scitbx::af::tiny<FloatType,2> q_range()
    {
      scitbx::af::tiny<FloatType,2> result;
      result[0] = min_q_;
      result[1] = max_q_;
      return(result);
    }

    scitbx::af::tiny<FloatType,3> binning_setup()
    {
      scitbx::af::tiny<FloatType,3> result;
      result[0]=start_q_;
      result[1]=end_q_;
      result[2]=q_step_;
      return(result);
    }

    private:
      int nx_, ny_, n_bins_;
      FloatType scalex_, scaley_, orgx_, orgy_, distance_, wavelength_, min_q_, max_q_, q_step_,start_q_, end_q_ ;
      scitbx::af::shared< FloatType > q_;
      scitbx::af::shared< FloatType > q_low_bin_values_;
      scitbx::af::shared< FloatType > q_mean_bin_values_;
      scitbx::af::shared< std::vector<int>  > q_low_bin_indices_;
      std::vector<int>  before_indices_, after_indices_;
  };




  template <typename FloatType>
  class image_mask
  {
    public:
      image_mask( int const& nx,
                  int const& ny )
      {
        nx_=nx;
        ny_=ny;

        SCITBX_ASSERT( nx > 0 );
        SCITBX_ASSERT( ny > 0 );

        // making the initial empty mask
        for (int ii=0;ii<nx_;ii++){ // looping over x
           for (int jj=0; jj<ny_; jj++){ // loop over y coordinate
             mask_.push_back( 0 );
           } //end looping over y
        } // end looping over x

      }

      void add_polygon(scitbx::af::const_ref<int> const& corners_x,
                       scitbx::af::const_ref<int> const& corners_y,
                       int const& inside_x, int const& inside_y)
      {
        inside_x_ = inside_x;
        inside_y_ = inside_y;
        SCITBX_ASSERT( corners_x.size() == corners_y.size() );
        SCITBX_ASSERT( corners_x.size() >= 3);

        for (int ii=0; ii<corners_x.size(); ii++){
          corners_x_.push_back( corners_x[ii] );
          corners_y_.push_back( corners_y[ii] );
        }
        fill_connecting_lines();
      }

      void add_circle( int center_x, int center_y, int radius )
      {

        radius=radius*radius;
        int tmp,count=0;
        for (int ii=0;ii<nx_;ii++){
          for (int jj=0;jj<ny_;jj++){
            tmp = (center_x-ii)*(center_x-ii) + (center_y-jj)*(center_y-jj);
            if (tmp <= radius){
              mask_[ count ] = 1;
            }
            count += 1;
          }
        }


      }


      void add( scitbx::af::shared< int > const& external_mask )
      {
        for (int ii=0;ii<external_mask.size();ii++){
          if (external_mask[ii]==1){
            mask_[ii]=1;
          }
        }

      }


      scitbx::af::shared< int > mask()
      {
        return(mask_);
      }

      void detect_dead_pixels( scitbx::af::shared< int > const& sample_image )
      {
        // a dead pixel is a pixel with value 0
        for (int ii=0;ii<sample_image.size();ii++){
          if (sample_image[ii]==0){
            mask_[ii]=1;
          }
        }
        // please dilate once to get rid of unreliable border pixels
        int n_changes;
        n_changes =mask_dilate();


      }

      int mask_dilate()
      {
         // loop over al points. when a mask point is encountered, set the neighbouring points to 9 if not 1 allready
         int count=0, a,b,c,d, n_changes=0;
         for (int ii=0;ii<nx_;ii++){
           for(int jj=0;jj<ny_;jj++){
             if (mask_[count]==1){
               a = count;
               b = count;
               c = count;
               d = count;
               if (ii>0){
                 a = x_y_to_count(ii-1,jj);
               }
               if (ii<nx_-1){
                 b = x_y_to_count(ii+1,jj);
               }
               if (jj>0){
                 c = x_y_to_count(ii,jj-1);
               }
               if (jj<ny_-1){
                 d = x_y_to_count(ii,jj+1);
               }
               if (mask_[a]==0){ mask_[a]=9; n_changes++; }
               if (mask_[b]==0){ mask_[b]=9; n_changes++; }
               if (mask_[c]==0){ mask_[c]=9; n_changes++; }
               if (mask_[d]==0){ mask_[d]=9; n_changes++; }
               count += 1; 
             }
           }
         }
         // all done with changes. replace the 9's with 1's
         mask_replace(9,1); 
      }

      void mask_invert()
      {
        mask_replace(0,9);
        mask_replace(1,0);
        mask_replace(9,1);
      }


    private:

      int x_y_to_count(int x, int y){
        int result;
        result = ny_*x + y;
        return ( result );
      }

      void fill_connecting_lines()
      {
        int count;
        int x1,x2,y1,y2;
        for (int kk=1;kk<corners_x_.size();kk++){
          x1 = corners_x_[kk-1];
          y1 = corners_y_[kk-1];
          x2 = corners_x_[kk];
          y2 = corners_y_[kk];
          set_line(x1,y1,x2,y2);
        }
        x1 = corners_x_[0];
        y1 = corners_y_[0];
        x2 = corners_x_[corners_x_.size()-1];
        y2 = corners_y_[corners_x_.size()-1];
        set_line(x1,y1,x2,y2);
        // -----------
        // set the growth point
        mask_[ x_y_to_count(inside_x_, inside_y_) ] = 9;
        int changes = 1;
        int a,b,c,d,m;
        while (changes > 0 ){
          changes = 0;
          // loop over all points;
          for (int ii=0;ii<nx_;ii++){
            for (int jj=0;jj<ny_;jj++){
               m = x_y_to_count(ii,jj);

               if (mask_[m]==9){
                 a = m;
                 b = m;
                 c = m;
                 d = m;


                 if (ii>0){
                   a = x_y_to_count(ii-1,jj);
                 }
                 if (ii<nx_-1){
                   b = x_y_to_count(ii+1,jj);
                 }
                 if (jj>0){
                   c = x_y_to_count(ii,jj-1);
                 }
                 if (jj<ny_-1){
                   d = x_y_to_count(ii,jj+1);
                 }
                 // sett all top 9 unless it is 1
                 if (mask_[a]==0){ mask_[a]=9; changes+=1;}
                 if (mask_[b]==0){ mask_[b]=9; changes+=1;}
                 if (mask_[c]==0){ mask_[c]=9; changes+=1;}
                 if (mask_[d]==0){ mask_[d]=9; changes+=1;}



               } /// was 9

            }
          }
        }
        mask_replace(0,1);
        mask_replace(9,0);
      }


      void mask_replace(int now, int later)
      {
         int kk;
         for (int ii=0;ii<nx_;ii++){
           for (int jj=0;jj<ny_;jj++){
             kk = x_y_to_count( ii, jj);
             if ( mask_[kk]==now){
               mask_[kk] = later;
             }
           }
         }
      }

      #define SWAP(x,y) (x ^= y ^= x^= y )
      void set_line(int x0, int y0, int x1, int y1)
      {
        bool order, steep;
        steep = std::abs(y1 - y0) > std::abs(x1 - x0);
        if (steep){ // we don't like steep lines
          SWAP(x0,y0);
          SWAP(x1,y1);
        }
        order = x0>x1;
        if (order){ // make sure that x0 < x1
          SWAP(x0,x1);
          SWAP(y0,y1);
        }
        int x,y;
        FloatType ft_x, ft_y,rc=0;
        if (x0 != x1){
          rc = static_cast<FloatType>(y0-y1)/static_cast<FloatType>(x0-x1);
          for (int ii=0;ii<=std::abs(x0-x1);ii++){
            x = x0 + ii;
            y = static_cast<int>( std::floor( y0 + ii*rc +0.5 ) );
            if (steep){
              mask_[ x_y_to_count(y,x) ]=1;
            } else {
              mask_[ x_y_to_count(x,y) ]=1;
            }
          }
        } else {
          x = x0;
          for (int ii=0;ii<std::abs(y0-y1);ii++){
            y = std::min(y0,y1)+ii;
            if (steep){
              mask_[ x_y_to_count(y,x) ]=1;
            } else {
              mask_[ x_y_to_count(x,y) ]=1;
            }
          }
        }
      }




     protected:
      scitbx::af::shared<int> mask_;
      std::vector<FloatType> corners_x_, corners_y_;
      int nx_, ny_, inside_x_, inside_y_;

  };

  template <typename FloatType>
  class radial_integrator
  {
    public:
      radial_integrator(  perpendicular_detector<FloatType> const& detector,   // for now, we only have a perpendicular detector, this should be a generic detector though in the future
                          scitbx::af::const_ref<int> const& mask,              // This mask defines the area of the detector that has no signal nor 'dark' info
                          scitbx::af::const_ref<int> const& shadow )           // This mask define the area of the detector that has 'dark' info
      :
      detector_( detector )
      {
         /*
             For clarity, lets define
             signal   pixel:  mask = 0, shadow = 0
             ignored  pixel:  mask = 1
             dark     pixel:  mask = 0, shadow = 1
         */
         SCITBX_ASSERT( (mask.size() == shadow.size()) || (shadow.size()==0) );
         for (int ii=0;ii<mask.size();ii++){
           mask_.push_back( mask[ii] );
           if (shadow.size()==0){
             shadow_.push_back( 0 );
           } else {
             shadow_.push_back( shadow[ii] );
           }
         }

         // setup the bins and indices
         bins_ = detector_.q_low_bin_values();
         indices_ = detector_.q_bin_indices();
         // init the moments array
         n_moments_ = 4;
         for (int ii=0;ii<bins_.size();ii++){
           std::vector< FloatType > tmp;
           for (int jj=0;jj<=n_moments_;jj++){
             tmp.push_back( 0.0 );
           }
           moments_.push_back( tmp );
         }

      }

      void integrate(scitbx::af::const_ref< int > const& data, bool shadow_only)
      {
         // first, set all bins to zero
         for (int ii_bin=0;ii_bin<bins_.size();ii_bin++){
           for (int kk=0;kk<moments_[ii_bin].size();kk++){
             moments_[ii_bin][kk]=0;
           }
         }


         FloatType tmp1, tmp2;
         // we have the indices, now we can perform the integration
         int index_in_image;
         bool go;
         for (int ii_bin=0;ii_bin<bins_.size();ii_bin++){
           for (int this_index=0;this_index<indices_[ ii_bin ].size();this_index++){
             // we have to make sure that we only add data (or shadow)
             index_in_image = indices_[ii_bin][this_index];

             go = false;
             if ( mask_[ index_in_image ]==0){ // this is not a masked pixel
               if (shadow_only){ // we only are interested in the shadow
                 if (shadow_[ index_in_image ]==1){
                   go = true;
                 }
               } else { // please ignore the masked area
                 if (shadow_[ index_in_image ]==0){
                   go = true;
                 }
               }
             }

             if (go){
               moments_[ ii_bin ][0] += 1; // count the contributors
               tmp1 = static_cast<FloatType>( data[ indices_[ii_bin][this_index] ] );
               tmp2 = tmp1;
               for (int kk = 1;kk<n_moments_+1; kk++){
                 moments_[ii_bin][kk]+=tmp2;
                 tmp2 = tmp1*tmp2;
               }
             }

           }
         }
         // do the averaging please
         for (int ii_bin=0;ii_bin<bins_.size();ii_bin++){
           for (int kk = 1;kk<n_moments_+1; kk++){
             moments_[ii_bin][kk]/=(moments_[ii_bin][0]+1e-13); // add small number to avoid div by zero
           }
         }
      }

      scitbx::af::shared<FloatType> mean_intensity()
      {
         scitbx::af::shared<FloatType> result;
         for (int ii=0;ii<moments_.size();ii++){
           result.push_back( moments_[ii][1] );
         }
         return(result);
      }

      scitbx::af::shared<FloatType> variance_intensity()
      {
         scitbx::af::shared<FloatType> result;
         for (int ii=0;ii<moments_.size();ii++){
           result.push_back( (moments_[ii][2]-moments_[ii][1]*moments_[ii][1])/(moments_[ii][0]+1e-13) );
         }
         return(result);
      }




    protected:
      perpendicular_detector<FloatType> detector_;
      scitbx::af::shared<int> mask_, shadow_;
      scitbx::af::shared< std::vector<FloatType> > moments_;
      scitbx::af::shared< std::vector<int> > indices_;
      scitbx::af::shared< FloatType > bins_;
      int n_moments_;
  };



}}  // namespace sastbx::data_reduction
#endif // SASTBX_DATA_REDUCTION_H
