//! Peter Zwart April 05, 2005
#ifndef SASTBX_INTENSITY_MODEL_H
#define SASTBX_INTENSITY_MODEL_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <map>

#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <string>
#include <iomanip>
#include <sastbx/intensity/scatlib.h>

namespace sastbx { namespace intensity {

  template <typename FloatType>
  class model
  {
    public:
      model(  scitbx::af::const_ref< scitbx::vec3<FloatType> > const& xyz, // coordinates
                     scitbx::af::const_ref< FloatType > const& radius,
                     scitbx::af::const_ref< FloatType > const& b_value,
                     scitbx::af::const_ref< FloatType > const& occupancy,
                     scitbx::af::const_ref< std::string > const& atom_types,
                     scattering_library< FloatType > const& scat_lib,
                     bool b_factor_on
                  ):
      scat_lib_(scat_lib), B_factor_(b_factor_on),max_r_(0.0)
      {
        SCITBX_ASSERT( xyz.size() > 0 );
        SCITBX_ASSERT( xyz.size() == b_value.size() );
        SCITBX_ASSERT( xyz.size() == occupancy.size() );

        n_ = xyz.size();
// Initializing local variables
        for(int i=0; i < n_; i++)
        {
         xyz_.push_back(xyz[i]);
         radius_.push_back(radius[i]);
         b_.push_back(b_value[i]);
         occ_.push_back(occupancy[i]);
         atom_types_.push_back(atom_types[i]);
        } // done with initialization

// (1) Re-Center the coordinates & Convert to Polar coordinates
        center_=Find_center();
        for(int i=0;i<n_;i++)
        {
                xyz_[i]-=center_;
                rtp_.push_back(Polar(i));
                if(rtp_[i][0] > max_r_) max_r_=rtp_[i][0];
        }
// (2) get scattering factors
        for(int i=0; i< n_; i++)
        {
        scattering_factor_indx_.push_back( scat_lib_.get_sf_indx(atom_types_[i]) );
        }
        time_t t;
        //cout<<time(&t)<<endl;
        Load_SF();
        //cout<<time(&t)<<endl;
      }//end of model constructor

      int size() {return n_;}

      FloatType distance(int i,int j){ return (xyz_[i]-xyz_[j]).length();}

      FloatType get_max_radius() {return max_r_;}
      scitbx::vec3<FloatType> get_xyz(int i) { return xyz_[i];}
      scitbx::vec3<FloatType> get_rtp(int i) { return rtp_[i];}
      scitbx::vec3<FloatType> get_center() { return center_;}

      scitbx::vec3<FloatType> Find_center()
      {
        scitbx::vec3<FloatType> center(0.0, 0.0, 0.0);
        for(int i=0;i<n_;i++)
        {
                center+=xyz_[i];
        }
        return center/n_;
      }

      scitbx::vec3<FloatType> Cart(scitbx::vec3<FloatType> polar)
      {
        scitbx::vec3<FloatType> this_xyz;
        this_xyz[0] = polar[0]* sin(polar[1]) * cos(polar[2]);
        this_xyz[1] = polar[0]* sin(polar[1]) * sin(polar[2]);
        this_xyz[2] = polar[0]* cos(polar[1]);
        return this_xyz;
      }

      scitbx::vec3<FloatType> Polar(int i)
      {
        scitbx::vec3<FloatType> plr(0,0,0);
        FloatType r,t,p;
        FloatType x,y,z;
        x=xyz_[i][0];
        y=xyz_[i][1];
        z=xyz_[i][2];
        r=xyz_[i].length();
        if(r==0.0) {
        return plr;
        }

        t=acos(z/r);
        p=atan2(y,x);
        plr[0]=r;
        plr[1]=t;
        plr[2]=p;
        return plr;
      }

      void Load_SF()
      {
        int n_q = scat_lib_.n_q();
        for(int index=0;index<n_q;index++)
        {
        scitbx::af::shared< FloatType > sf_i_array;
        scitbx::af::shared< FloatType > dummy_sf_i_array;
        sf_i_array.reserve(n_);
        dummy_sf_i_array.reserve(n_);

         for(int i=0;i<n_;i++)
         {
            if(B_factor_)
            {   sf_i_array[i]=scat_lib_.get_sf_with_b_and_occ(scattering_factor_indx_[i],index, b_[i],occ_[i]);
                dummy_sf_i_array[i]=scat_lib_.get_dummy_sf(scattering_factor_indx_[i],index);
            }
            else
            {   sf_i_array[i]=scat_lib_.get_sf(scattering_factor_indx_[i],index);
                dummy_sf_i_array[i]=scat_lib_.get_dummy_sf(scattering_factor_indx_[i],index);
            }
         }
         sf_array_.push_back(sf_i_array);
         dummy_sf_array_.push_back(dummy_sf_i_array);
        }
        return;
      }

    FloatType get_radius(int i)
    {
        return radius_[i];
    }
    public:
      scitbx::af::shared< scitbx::af::shared<FloatType> > sf_array_;
      scitbx::af::shared< scitbx::af::shared<FloatType> > dummy_sf_array_;

    private:
      scitbx::af::shared< scitbx::vec3<FloatType> > xyz_;
      scitbx::af::shared< scitbx::vec3<FloatType> > rtp_;
      scitbx::vec3<FloatType> center_;
      scitbx::af::shared< FloatType > radius_;
      scitbx::af::shared< FloatType > b_;
      scitbx::af::shared< FloatType >  occ_;
      scitbx::af::shared< std::string > atom_types_;//same as scat_types_ in the lib, stored for each atom
      scattering_library< FloatType > scat_lib_;
      scitbx::af::shared< int > scattering_factor_indx_;
      int n_;
      FloatType max_r_;
      bool B_factor_;
  };


}}  // namespace sastbx::intensity
#endif // SASTBX_INTENSITY_MODEL_H
