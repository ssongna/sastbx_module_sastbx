//! Peter Zwart April 05, 2005
#ifndef SASTBX_INTENSITY_SCATLIB_H
#define SASTBX_INTENSITY_SCATLIB_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <map>

#include <scitbx/vec3.h>
#include <string>
#include <iomanip>


namespace sastbx { namespace intensity {

  template <typename FloatType>
  class scattering_library
  {
    public:
      scattering_library( scitbx::af::const_ref< FloatType > const& q_values )
      {
        SCITBX_ASSERT( q_values.size() > 0 );
        for (int ii=0;ii<q_values.size();ii++){
          q_.push_back( q_values[ii] );
        }
        n_ = q_.size();
      }

    void load_scattering_info( std::string const& scat_type, scitbx::af::const_ref< FloatType > sf, scitbx::af::const_ref< FloatType > dummy_sf )
    {
      scat_types_.push_back( scat_type );
      scitbx::af::shared< FloatType > tmp_sf;
      scitbx::af::shared< FloatType > tmp_dummy_sf;
      SCITBX_ASSERT( sf.size() == q_.size() );
      for (int ii=0;ii<q_.size();ii++){
        tmp_sf.push_back( sf[ii] );
        tmp_dummy_sf.push_back( dummy_sf[ii] );
      }
      sfs_.push_back( tmp_sf );
      dummy_sfs_.push_back( tmp_dummy_sf );
    }


    FloatType q(int ii)
    {
      SCITBX_ASSERT( ii >= 0);
      SCITBX_ASSERT( ii < n_ );
      return(q_[ii]);
    }

    scitbx::af::shared< FloatType > q_range()
    {
      return ( q_ );
    }

    int n_q()
    {
      return ( q_.size() );
    }

    scitbx::af::shared <FloatType> get_sfs(std::string scat_type)
    {
        bool found_it=false;
        for (int i=0;i<scat_types_.size();i++)
        { if(scat_type == scat_types_[i]){
                found_it=true;
                return sfs_[i];
          }
        }
       SCITBX_ASSERT( found_it );//not found
    }

    int get_sf_indx(std::string scat_type)
    {
        bool found_it=false;
        for(int i=0; i<scat_types_.size();i++)
        {
         if(scat_type == scat_types_[i]) {
                found_it=true;
                return i;
          }
        }
        SCITBX_ASSERT( found_it ); //not found
    }

    FloatType get_sf( int ii, int qq)
    {
        SCITBX_ASSERT( ii < sfs_.size() );
        SCITBX_ASSERT( ii >= 0 );
        SCITBX_ASSERT( qq >= 0 );
        SCITBX_ASSERT( qq < q_.size() );
        return sfs_[ii][qq];
    }

    FloatType get_dummy_sf( int ii, int qq)
    {
        SCITBX_ASSERT( ii < sfs_.size() );
        SCITBX_ASSERT( ii >= 0 );
        SCITBX_ASSERT( qq >= 0 );
        SCITBX_ASSERT( qq < q_.size() );
        return dummy_sfs_[ii][qq];
    }

    FloatType get_sf_with_b_and_occ( int ii, int qq, FloatType b, FloatType occ)
    {
      FloatType result=get_sf(ii,qq);
      FloatType stol = q_[qq]/(scitbx::constants::pi*4.0);
      result = result*std::exp( -b*stol*stol )*occ;
      return result;
    }

    private:
      scitbx::af::shared< FloatType > q_;
      scitbx::af::shared< scitbx::af::shared<FloatType> > sfs_;
      scitbx::af::shared< scitbx::af::shared<FloatType> > dummy_sfs_;
      scitbx::af::shared< std::string > scat_types_;
      int n_;
  };

}}  // namespace sastbx::intensity
#endif // SASTBX_INTENSITY_SCATLIB_H
