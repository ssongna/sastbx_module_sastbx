#ifndef SASTBX_BEADS_H
#define SASTBX_BEADS_H
#endif

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>

#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <string>
#include <iomanip>
#include <tntbx/include/tnt.h>

namespace sastbx { namespace cgmd {

   template<typename FloatType>
   class beads
   {
     public:
       beads(scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& xyz,
               scitbx::af::const_ref <int> const& bead_indx,
               FloatType const& cutoff,
               FloatType const& sf): 
      cutoff_(cutoff), n_bead_(bead_indx.size() ),
      natom_( xyz.size() ),
      dist0_(n_bead_,n_bead_, 0.0), 
      dist1_(n_bead_,n_bead_, 0.0), 
      distE_(n_bead_,n_bead_, 0.0), 
      bead_xyz_(n_bead_,scitbx::vec3<FloatType>(0.0) ), 
      dxyz_(n_bead_,scitbx::vec3<FloatType>(0.0) ), 
      correct_xyz_(n_bead_,scitbx::vec3<FloatType>(0.0) ), 
      cws_(n_bead_,0.0),
      all_dxyz_(natom_,scitbx::vec3<FloatType>(0,0,0)),
      sf_(n_bead_,n_bead_, 0.0), sf_13_(sf), sf_14_(sf/2.0),
      gamma_(0.9), kamma_(0.05)
      {
         for(int i=0;i<n_bead_;i++)
         {
           bead_indx_.push_back(bead_indx[i]);
           bead_xyz_[i] = xyz[ bead_indx[i] ];
         }
         for(int i=0;i<natom_;i++)
           xyz_.push_back( xyz[i] );

         calc_dist();
         initialize_d0();
      }

      void update( scitbx::af::const_ref< scitbx::vec3< FloatType > > const& xyz )
      {
        for(int i=0;i<natom_;i++)
          xyz_[i] = xyz[i];
        calc_dist();
        update_de();
        return;
      }

      void initialize_d0()
      {
        for(int i=0;i<n_bead_;i++)
          for(int j=i+1;j<n_bead_;j++)
          {
             dist0_[i][j] =dist0_[j][i]= dist1_[i][j];
             distE_[i][j] =distE_[j][i]= dist1_[i][j];
             if( abs(i-j)<3 ) sf_[i][j]=sf_[j][i]=sf_13_;
             else if( abs(i-j)<4 ) sf_[i][j]=sf_[j][i]=sf_14_;
          }
      }

      void update_de()
      {
        FloatType fac_g(1-gamma_), fac_k(1-kamma_);
        for(int i=0;i<n_bead_;i++)
          for(int j=i+1;j<n_bead_;j++)
          {
             distE_[i][j] = fac_k*distE_[i][j]+kamma_*(gamma_*dist1_[i][j]+fac_g*dist0_[i][j]);
             distE_[j][i] = distE_[i][j];
          }

        return; 
      }

      void calc_dist( )
      {
         FloatType cutoff2(cutoff_*cutoff_);
         scitbx::vec3 <FloatType> vecDiff;
         FloatType distsqr, dist;
         FloatType scale_factor;
         dist1_ = 0.0;

         for(int i=0;i<n_bead_;i++)
         {
           for(int j=i+1;j<n_bead_;j++)
           {
             vecDiff = bead_xyz_[i]-bead_xyz_[j];
             distsqr = vecDiff * vecDiff;
             dist = sqrt( distsqr );
             if( distsqr == 0 ) { std::cout<<"distance between atoms is zero!("<<i<<","<<j<<")"<<std::endl; }
             if(distsqr <= cutoff2)
             {
               dist1_[i][j] = dist1_[j][i] = dist;
             }
           }
         }//end of for-loop n_bead_
      }

      FloatType restraint(scitbx::af::const_ref<scitbx::af::tiny<std::size_t,2> > indices,
                          scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& dxyz)
      {
        int np = indices.size();
        int i,j;
        FloatType d0, E(0.0), dist, sf;
        E = internal_restraint(dxyz);
        for(int n=0;n<np;n++)
        {
          i = indices[n][0];
          j = indices[n][1];
          d0 = distE_[i][j];
          sf = sf_[i][j];
          if( d0 == 0 ) continue;
          dist = (bead_xyz_[i]-bead_xyz_[j]+dxyz[i]-dxyz[j]).length();
          E += std::pow( (d0-dist), 2.0 )*sf;
        }
        std::cout<<E<<" "<<(n_bead_*3+np)<<std::endl;
        return E/static_cast<FloatType>(n_bead_*3+np);
      }

      FloatType internal_restraint( scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& dxyz ) 
      {
        int i,j;
        FloatType d12, d13, d14, E(0.0);
        for(i=0;i<n_bead_-3;i++)
        {
          d12 = (bead_xyz_[i]-bead_xyz_[i+1]+dxyz[i]-dxyz[i+1]).length();
          d13 = (bead_xyz_[i]-bead_xyz_[i+2]+dxyz[i]-dxyz[i+2]).length();
          d14 = (bead_xyz_[i]-bead_xyz_[i+3]+dxyz[i]-dxyz[i+3]).length();
          E += std::pow( (distE_[i][i+1]-d12),2.0)*sf_13_;
          E += std::pow( (distE_[i][i+2]-d13),2.0)*sf_13_;
          E += std::pow( (distE_[i][i+3]-d14),2.0)*sf_14_;
        }
        return E;
      }


// Relax the structure //
      scitbx::af::shared< scitbx::vec3< FloatType > > 
      relax( scitbx::af::const_ref<scitbx::af::tiny<std::size_t,2> > indices,
             scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& random_dxyz)
      {
        /*
        new_xyz = update_xyz()
        relax new_xyz to minimize internal_restraint to bearable level
        return new_dxyz for all-atom model update
        */
        scitbx::vec3<FloatType> d_vector(0.0);
        scitbx::vec3<FloatType> correction(0.0);
        FloatType d_new, d0, sf;
        Damp_vector( random_dxyz );
        int i,j,n, bi,bj;
        for(bi=0;bi<n_bead_;bi++)
          correct_xyz_[bi] = cws_[bi]=0.0;
        for(bi=0;bi<n_bead_-3;bi++)
        {
	  for(bj=bi+1;bj<bi+4;bj++)
          {
            d0 = dist0_[bi][bj];
            sf = sf_[bi][bj];
	    d_vector = (bead_xyz_[bi]-bead_xyz_[bj]+dxyz_[bi]-dxyz_[bj]);
            d_new = d_vector.length();
            correction = d_vector*( (d_new-d0)/d_new/2.0 );
            //std::cout<<d0<<" "<<d_new<<std::endl;
            correct_xyz_[bi] -= correction*sf;
            correct_xyz_[bj] += correction*sf;
            cws_[bi] += sf;
            cws_[bj] += sf;
          }
        } 

        int np = indices.size();
        for(n=0;n<np;n++)
        {
	  bi = indices[n][0];
          bj = indices[n][1];
          d0 = dist0_[bi][bj];
          if( d0 == 0 ) continue;
	  d_vector = (bead_xyz_[bi]-bead_xyz_[bj]+dxyz_[bi]-dxyz_[bj]);
          d_new = d_vector.length();
          correction = d_vector*( (d_new-d0)/d_new/2.0 );
          correct_xyz_[bi] += correction;
          correct_xyz_[bj] -= correction;
          cws_[bi] +=1.0;
          cws_[bj] +=1.0;
        }
        // add correction to dxyz //
        for(n=0;n<n_bead_;n++)
        {
          dxyz_[n] += ( correct_xyz_[n]/cws_[n] ) ;
        }

        return dxyz_;
      }

      scitbx::af::shared<scitbx::vec3< FloatType > >Project2All( 
        scitbx::af::const_ref< scitbx::vec3< FloatType > > dxyz)
      {
        for(int j=0;j<bead_indx_[1];j++)
          all_dxyz_[j] = dxyz[0];
        for(int i=1;i<n_bead_-1;i++)
          for(int j=bead_indx_[i];j<bead_indx_[i+1];j++)
            all_dxyz_[j] = dxyz[i];

        int bead_last=bead_indx_[n_bead_-1];
        for(int j=bead_last;j<natom_;j++)
          all_dxyz_[j] = dxyz[n_bead_-1];
        return all_dxyz_;
      }

      void Damp_vector(
        scitbx::af::const_ref< scitbx::vec3<FloatType> >const& dxyz )
      {
        int np = dxyz.size();
        scitbx::vec3< FloatType > center(0.0);
/*
        dxyz_[0] = (dxyz[0]+dxyz[1])/2.0;
        dxyz_[np-1] = (dxyz[np-1]+dxyz[np-2])/2.0;
        center = dxyz_[0]+dxyz_[np-1];
        for(int i=1;i<np-1;i++)
        {
          dxyz_[i] = (dxyz[i-1]+dxyz[i]+dxyz[i+1])/3.0;
          center += dxyz_[i];
        }
*/
        for(int i=0;i<np;i++)
        {
          dxyz_[i] = dxyz[i];
          center += dxyz_[i];
        }
        center /= static_cast<FloatType>(np);
        for(int i=0;i<np;i++)
          dxyz_[i] -= center;
        return;
      }

      scitbx::af::shared< FloatType > Damp_vector( scitbx::af::const_ref < FloatType > dx )
      {
        FloatType upper, lower;
        scitbx::af::shared< FloatType > new_dx( dx.size(),0 );
        int dimen = dx.size()-3, j;
        scitbx::vec3<FloatType> center(0.0);

        for(j=0;j<3;j++) //front and end averaged to the nearest nbr
        {
          new_dx[j] = (dx[j]+dx[j+3])/2.0;
          new_dx[j+dimen] = (dx[dimen+j]+dx[dimen+j-3])/2.0;
          center[j] += new_dx[j];
          center[j] += new_dx[dimen+j];
        }
        for(int i=3;i<dimen;i+=3)
        {
          for(j=0;j<3;j++)
          {
            new_dx[i+j] = (dx[i+j-3]+dx[i+j]+dx[i+j+3])/3.0;
            center[j] += new_dx[i+j];
          }
        }
        center /= static_cast<FloatType>( dimen/3+1 );

        for(int i=0;i<dimen+3;i+=3)
        {
          for(j=0;j<3;j++)
            new_dx[i+j] -= center[j];
        }

        return new_dx;
      }

      private:
        FloatType cutoff_,  sf_13_, sf_14_;
        int natom_, n_bead_;
        FloatType gamma_, kamma_;
        scitbx::af::shared <int> bead_indx_;
        scitbx::af::shared <scitbx::vec3< FloatType> > all_dxyz_;
        scitbx::af::shared <scitbx::vec3< FloatType> > xyz_;
        scitbx::af::shared <scitbx::vec3< FloatType> > bead_xyz_;
        scitbx::af::shared <scitbx::vec3< FloatType> > dxyz_;
        scitbx::af::shared <scitbx::vec3< FloatType> > correct_xyz_;
        scitbx::af::shared <FloatType> cws_; //correction_weight_sum
        TNT::Array2D <FloatType> dist0_;  // initial dist
        TNT::Array2D <FloatType> dist1_;  // instance dist
        TNT::Array2D <FloatType> distE_;  // Equilibrum dist
        TNT::Array2D <FloatType> sf_;
   }; //end of bead class


} //namepsace cgmd
} //namespace sastbx
