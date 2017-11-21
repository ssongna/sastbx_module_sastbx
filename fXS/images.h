#ifndef SASTBX_2dfXS_Image_H
#define SASTBX_2dfXS_Image_H

#include <cmath>
#include <iostream>
#include <scitbx/constants.h>
#include <scitbx/array_family/shared.h>
#include<scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/sort.h>

#include <scitbx/vec3.h>
#include <scitbx/vec2.h>

namespace array=scitbx::af;

namespace sastbx { namespace fXS {

  template <typename FloatType>
  class image_cartesian
  {
    public:
      image_cartesian() {}
      image_cartesian(
        int const& nx,
        int const& ny,
        array::const_ref< scitbx::vec3< FloatType > > image //(x_indx,y_indx,color)
      ): nx_(nx), ny_(ny)
      {
        int total_pix = image.size();
        FloatType r,t, xx,yy;
        center_x_ = (nx_-1.0)/2.0;
        center_y_ = (ny_-1.0)/2.0;
        radius_ = center_x_ ; //std::sqrt( center_x_*center_x_ + center_y_*center_y_ );
        for(int i=0; i<total_pix; i++)
        {
          xx = image[i][0]-center_x_;
          yy = image[i][1]-center_y_;
          r=std::sqrt( xx*xx + yy*yy );
          if( r<=radius_ )
          {
            t=std::atan2( yy,xx );
            if( t<0 ) t+=scitbx::constants::two_pi;
            rt_.push_back( scitbx::vec2<FloatType> ( r,t ) ) ;
            pixels_.push_back( image[i][2] );
          }
          sp_image_.push_back( image[i] );
        }
      }

      void calc_c2_array(int const& nr, int const& nt, bool further_reduct, bool smoothen)
      {
        setup_polar_bin( nr, nt, further_reduct );
        int nx(nr*2+1), ny(nr+1);
        setup_carte_bin( nx, ny );
        FloatType dr = ( radius_/(nr) );
        FloatType dt = ( scitbx::constants::two_pi / nt );
        int i,j, nt_calc;
        int x,y;
        FloatType r, t, this_c2, sasI_sqr;
        nt_calc = nt/2;  // Only a quarter of plain needs to be evaluated


        for(i=3;i<nr;i++) //skips the small q's
        {
          r=dr*i;
          for(j=0;j<=nt_calc;j++)  //the lower half can be obtained by mirror symmetry
          {
            t=dt*j;
            x=int( r*std::cos(t) +nr );
            y=int( r*std::sin(t) );
            this_c2 = compute_c2( i,j,nt, bin_image_[i]);
            xy_bin_image_[x][y] += this_c2;
            xy_bin_count_[x][y] += 1.0;
          }
        }

        int n_tot( 2*nr+1 );
        array::c_grid<2> value_grid( n_tot, n_tot );
        array::versa<  FloatType, array::c_grid<2> > values(value_grid,0);
        for(i=0;i<nx;i++)
        {
          for(j=0;j<ny;j++)
          {
            this_c2=0.0;
            if(xy_bin_image_[i][j] !=0 )
            {
              this_c2 = xy_bin_image_[i][j]/xy_bin_count_[i][j];
            }
            values(i,nr+j) = this_c2;
            if( j>0 )
              values(i,nr-j) = this_c2;
          }
        }
        int radius = 1;
        if(smoothen)
          values = median_filter( radius, values.const_ref() );

        c2_image_.clear();
        for(int i=0;i<n_tot;i++)
          for(int j=0;j<n_tot;j++)
            c2_image_.push_back( scitbx::vec3< FloatType> (i,j,values(i,j) ) );

      }
        

/*            c2_image_.push_back( scitbx::vec3< FloatType >(i,nr+j,this_c2) );

            if( i>0 && j>0 )
            {
              c2_image_.push_back( scitbx::vec3< FloatType >(nr-i,nr-j,this_c2) );
              c2_image_.push_back( scitbx::vec3< FloatType >(nr+i,nr-j,this_c2) );
              c2_image_.push_back( scitbx::vec3< FloatType >(nr-i,nr+j,this_c2) );
            }
            else
            {
            if( j>0)
              c2_image_.push_back( scitbx::vec3< FloatType >(i,nr-j,this_c2) );
//            if( i>0)
//              c2_image_.push_back( scitbx::vec3< FloatType >(nr-i,nr+j,this_c2) );
//            }
          }
        }
*/

      scitbx::af::versa<FloatType, scitbx::af::c_grid<2>  > median_filter(int const& radius, array::const_ref<FloatType, array::c_grid<2> > values)
      {
         int n_tot=int(std::sqrt( values.size() ) );
         int median_point = static_cast<int>((radius*2+1)*(radius*2+1)/2.0+0.5); // technically not correct, but sufficient enough

         scitbx::af::c_grid<2> value_grid( n_tot, n_tot );
         scitbx::af::versa<  FloatType, scitbx::af::c_grid<2> > new_value(value_grid,0);
         // now we have to walk over the whole image
         for (int xx=0+radius;xx<n_tot-radius;xx++){
           for (int yy=0+radius;yy<n_tot-radius;yy++){
             scitbx::af::shared<FloatType> tmp_vals;
             for (int ii=-radius;ii<radius+1;ii++){
               for (int jj=-radius;jj<radius+1;jj++){
                  tmp_vals.push_back( values(xx+ii,yy+jj) );
               }
             }
             scitbx::af::shared<std::size_t> permut;
             permut = scitbx::af::sort_permutation( tmp_vals.const_ref(), true );
             FloatType median = tmp_vals[ permut[median_point] ];
             new_value( xx,yy ) = median;
           }
         }
         return(new_value);
      }


      array::shared< scitbx::vec3< FloatType > > get_c2_array()
      {
        if( c2_image_.size() == 0 )
          calc_c2_array(50,40,true,true); //calculate c2 pattern with default resolution
        return c2_image_;
      }

      FloatType compute_c2( int const& r_indx, int const& t_diff, int const& nt, array::shared<array::shared<int> > r_array )
      {
        FloatType this_c2(0.0);
        int nt_calc = nt/2;
        int nt0, nt1, t0, t1; //number of points in two theta_bin's and the indices
        int npairs(0);
        for(int i=0;i<nt_calc;i++)
        {
          nt0 = r_array[i].size();
          nt1 = r_array[i+t_diff].size();
          for( t0 = 0; t0<nt0; t0++ )
            for( t1 = 0; t1<nt1; t1++ )
            {
               this_c2 += pixels_[ r_array[i][t0] ]*pixels_[ r_array[i+t_diff][t1] ];
            }
          npairs += nt0*nt1;
/*
          nt0 = r_array[nt-i].size();
          nt1 = r_array[nt-i-t_diff].size();
          for( t0 = 0; t0<nt0; t0++ )
            for( t1 = 0; t1<nt1; t1++ )
            {
               this_c2 += pixels_[ r_array[nt-i][t0] ]*pixels_[ r_array[nt-i-t_diff][t1] ];
            }
          npairs += nt0*nt1;
*/
        }
        if( npairs == 0 ) return 0.0;
        return this_c2/static_cast< FloatType >(npairs);
      }

      void setup_polar_bin(int const& nr, int const& nt, bool further_reduct)
      {
        FloatType dr = ( radius_/(nr) );
        FloatType dt = ( scitbx::constants::two_pi / nt );
        FloatType Isas;
        int total_pix = rt_.size();
        int r_indx, t_indx, n_pixel_in_bin;
        int i,j;
        scitbx::vec2<FloatType> this_rt;

        bin_image_.clear();
        bin_r_.clear();
        bin_rcount_.clear();
        for(i=0;i<nr+1;i++)
        {
          array::shared< array::shared<int> > r_array;
          for(j=0;j<=nt;j++) {
            array::shared< int > rt_array;
            r_array.push_back( rt_array );
          }
          bin_image_.push_back( r_array );
          bin_r_.push_back( 0.0 );
          bin_rcount_.push_back( 0.0 );
        }

        for(i=0;i<total_pix; i++)
        {
          this_rt = rt_[i];
          r_indx=std::floor( this_rt[0]/dr + 0.5 );
          t_indx=std::floor( this_rt[1]/dt + 0.5 );
          if( pixels_[i] != 0 ) {
            bin_image_[r_indx][t_indx].push_back(i); 
            bin_r_[r_indx] += pixels_[i];
            bin_rcount_[r_indx] += 1;
          }
        }

        if( further_reduct )
        { 
        for(r_indx=0;r_indx<=nr;r_indx++)
        {
          if( bin_r_[r_indx] !=0 )
            bin_r_[r_indx] /= bin_rcount_[r_indx];
          Isas = bin_r_[r_indx];
          std::cout<<r_indx<<" "<<Isas<<" SAXS"<<std::endl;
          if( Isas > 0 ) 
          {
            for(t_indx=0;t_indx<=nt;t_indx++)
            {
              n_pixel_in_bin = bin_image_[r_indx][t_indx].size();
              for( i=0;i<n_pixel_in_bin;i++)
              {
                pixels_[bin_image_[r_indx][t_indx][i] ] -= Isas; 
                pixels_[bin_image_[r_indx][t_indx][i] ] /= Isas; 
              }
            }
          }
        }
        for(i=0;i<total_pix;i++)
          std::cout<<rt_[i][0]*std::cos(rt_[i][1])<<" "<<rt_[i][0]*std::sin(rt_[i][1])<<" "<<pixels_[i]<<std::endl;
        normalize_sp( dr );
        } //end of further deduction (subtract radial average and normalize)

        return;
     }

     array::shared< scitbx::vec3< FloatType > > get_normalized_sp() 
     {
       return sp_image_;
     } 
  
     void normalize_sp( FloatType dr)
     {
        int n_pixel = sp_image_.size();
        int r_indx;
        FloatType xx,yy,r;
        for(int i=0;i<n_pixel;i++)
        {
          xx = sp_image_[i][0]-center_x_;
          yy = sp_image_[i][1]-center_y_;
          r=std::sqrt( xx*xx + yy*yy );
          if( r<=radius_ )
          {
            r_indx = std::floor(r/dr);
            if( bin_r_[r_indx] == 0 ) sp_image_[i][2] = 0.0;
            else  
            { sp_image_[i][2] -= bin_r_[r_indx];
              sp_image_[i][2] /= bin_r_[r_indx];
            }
          }
          else sp_image_[i][2] = 0.0;
        }
        return;
      }

      void setup_carte_bin( int const& nx, int const& ny )
      {
        xy_bin_image_.clear();
        xy_bin_count_.clear();
        for(int i=0;i<nx;i++)
        {
          array::shared< FloatType > t_array( ny, 0.0);
          array::shared< FloatType > t_array_count( ny, 0);
          xy_bin_image_.push_back( t_array );
          xy_bin_count_.push_back( t_array_count );
        }
        return;
      }

    private:
      int nx_, ny_;
      FloatType center_x_, center_y_, radius_;
      array::shared< FloatType > pixels_;
      array::shared< FloatType > bin_r_;
      array::shared< FloatType > bin_rcount_;
      array::shared< scitbx::vec2< FloatType > > rt_;
      array::shared< array::shared< array::shared<int> > > bin_image_;
      array::shared< array::shared< FloatType > > xy_bin_image_;
      array::shared< array::shared< FloatType > > xy_bin_count_;
      array::shared< scitbx::vec3< FloatType > > sp_image_;
      array::shared< scitbx::vec3< FloatType > > c2_image_;
  };

  template <typename FloatType>
  class image_polar
  {
    image_polar() {}
  };

}} //namespace sastbx::fXS


#endif
