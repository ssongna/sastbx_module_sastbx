#ifndef SASTBX_ELASTIC_H
#define SASTBX_ELASTIC_H
#endif

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <scitbx/histogram.h>
#include <map>

#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <string>
#include <iomanip>
#include <tntbx/include/tnt.h>
#include <tntbx/include/jama_eig.h>

namespace sastbx { namespace refine {

   template<typename FloatType>
   class elastic
   {
        public:
          elastic(scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& xyz,
                  scitbx::af::const_ref <int> const& ca_indx) //CA based elastic network
          {
           FloatType cutoff2(cutoff_*cutoff_);
           natom_ = xyz.size();
           n_ca_  = ca_indx.size();
           mean_disp_=1.0/sqrt(static_cast<FloatType> (n_ca_*3));
           dist_array.reserve(n_ca_*(n_ca_-1)/2);
           for(int i=0;i<n_ca_*(n_ca_-1)/2;i++)
                dist_array.push_back(0);

           for(int i=0;i<n_ca_;i++)
                ca_indx_.push_back(ca_indx[i]);
          }
//Elastic Network Model Plus Bead Restraint functionality //
//
          elastic(scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& xyz,
                  scitbx::af::const_ref <int> const& ca_indx, //CA based elastic network
                  FloatType const& cutoff,
                  FloatType const& sf): cutoff_(cutoff), sf_(sf),
          natom_( xyz.size() ),
          n_ca_(ca_indx.size()),
          Hessian(n_ca_*3, n_ca_*3, 0.0),
          dist0_(n_ca_, n_ca_, 0.0),
          this_dxyz_(n_ca_,scitbx::vec3<FloatType>(0.0) ),
          correct_xyz_(n_ca_,scitbx::vec3<FloatType>(0.0) ),
          cws_(n_ca_,0.0)
          {
           FloatType cutoff2(cutoff_*cutoff_);
           mean_disp_=1.0/sqrt(static_cast<FloatType> (n_ca_*3));
           dist_array.reserve(n_ca_*(n_ca_-1)/2);
           for(int i=0;i<n_ca_*(n_ca_-1)/2;i++)
             dist_array.push_back(0);

           for(int i=0;i<n_ca_;i++)
           {
             ca_indx_.push_back(ca_indx[i]);
             ca_xyz_.push_back( xyz[ca_indx[i] ] );
           }

           int hi,hj,indx_i,indx_j;
           scitbx::vec3 <FloatType> vecDiff;
           FloatType distsqr, Hessij;
           FloatType scale_factor;

           for(int i=0;i<n_ca_;i++)
           { hi = i*3;
             for(int j=i+1;j<n_ca_;j++)
             {  hj=j*3;
                indx_i=ca_indx[i];
                indx_j=ca_indx[j];
                vecDiff = xyz[indx_i]-xyz[indx_j];
                distsqr = vecDiff * vecDiff;
                if( distsqr == 0 ) { std::cout<<"distance between atoms is zero!("<<indx_i<<","<<indx_j<<")"<<std::endl; }
                if(j-i < 4) { scale_factor = sf_;}
                else { scale_factor = 1;}
                if(distsqr <= cutoff2)
                {
                  dist0_[i][j] = dist0_[j][i] = sqrt( distsqr );  //distance matrix for restraint calculation
                  distsqr *= scale_factor;
                  for(int ii=0;ii<3;ii++)
                    for(int jj=ii;jj<3;jj++)
                    {
                        Hessij = vecDiff[ii]*vecDiff[jj]/distsqr; //*scale_factor;
                        if(ii == jj) Hessij /= 2.0;
                        Hessian[hi+ii][hi+jj] += Hessij;
                        Hessian[hi+ii][hj+jj] -= Hessij;
                        Hessian[hi+jj][hi+ii] += Hessij;
                        Hessian[hi+jj][hj+ii] -= Hessij;
                        Hessian[hj+ii][hj+jj] += Hessij;
                        Hessian[hj+jj][hj+ii] += Hessij;
                        Hessian[hj+ii][hi+jj] -= Hessij;
                        Hessian[hj+jj][hi+ii] -= Hessij;
                    }
                }
             }
           }//end of for-loop n_ca_
//        std::cout<<ca_indx.size()<<" "<<xyz.size()<<std::endl;
//        std::cout<<Hessian;
          JAMA::Eigenvalue <FloatType> eig(Hessian);
          eig.getV(eigenvec);  // ! row oriented, but each column is one eigenvector
          eig.getRealEigenvalues(eigenval);
        //std::cout << eigenval;
//        std::cout<<eigenvec;
          }

        int nmode(FloatType percentage, int n_max_mode)
        {
          FloatType total(0.0), partial(0.0);
          int n_val( eigenval.dim() ), i;
          for(i=6;i<n_val;i++)
            total += 1.0/eigenval[i];
          total *= percentage;
          for(i=6;i<n_max_mode;i++)
          {
            partial += 1.0/eigenval[i];
            if( partial >= total) break;
          }
          return i-6;
        }

        scitbx::af::shared<FloatType> eigenvalues()
        {
          scitbx::af::shared<FloatType> amp;
          int n_val( eigenval.dim() );
          for(int i=6;i<n_val;i++)
            amp.push_back(1.0/(eigenval[i]+1e-18));
          return amp;
        }

      FloatType restraint(scitbx::af::const_ref<scitbx::af::tiny<std::size_t,2> > indices,
                          scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& dxyz)
      {
        SCITBX_ASSERT( dxyz.size() == n_ca_ );
        int np = indices.size();
        int i,j;
        FloatType d0, E(0.0), internal_E, dist;
        internal_E = internal_restraint(dxyz);
        for(int n=0;n<np;n++)
        {
          i = indices[n][0];
          j = indices[n][1];
          d0 = dist0_[i][j];
          if( d0 == 0 ) continue;
          dist = (ca_xyz_[i]-ca_xyz_[j]+dxyz[i]-dxyz[j]).length();
          E += std::pow( (d0-dist), 2.0 );
        }
        E += internal_E*n_ca_*3;
        return E/static_cast<FloatType>(n_ca_*3+np);
      }

      FloatType internal_restraint( scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& dxyz )
      {
        int i,j;
        FloatType d12, d13, d14, E(0.0);
        for(i=0;i<n_ca_-3;i++)
        {
          d12 = (ca_xyz_[i]-ca_xyz_[i+1]+dxyz[i]-dxyz[i+1]).length();
          d13 = (ca_xyz_[i]-ca_xyz_[i+2]+dxyz[i]-dxyz[i+2]).length();
          d14 = (ca_xyz_[i]-ca_xyz_[i+3]+dxyz[i]-dxyz[i+3]).length();
          E += std::pow( (dist0_[i][i+1]-d12),2.0);
          E += std::pow( (dist0_[i][i+2]-d13),2.0);
          E += std::pow( (dist0_[i][i+3]-d14),2.0);
        }
        return E/static_cast<FloatType>(n_ca_*3);
      }



// Relax the structure //
      scitbx::af::shared< scitbx::vec3< FloatType > >
      relax( scitbx::af::const_ref<scitbx::af::tiny<std::size_t,2> > indices,
             scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& dxyz)
      {
        /*
        new_xyz = update_xyz()
        relax new_xyz to minimize internal_restraint to bearable level
        return new_dxyz for all-atom model update
        */
        scitbx::vec3<FloatType> d_vector(0.0);
        scitbx::vec3<FloatType> correction(0.0);
        FloatType d_new, d0, sf;

        int i,j,n, bi,bj;
        for(bi=0;bi<n_ca_;bi++)
          correct_xyz_[bi] = cws_[bi]=0.0;
        for(bi=0;bi<n_ca_-3;bi++)
        {
          for(bj=bi+1;bj<bi+4;bj++)
          {
            d0 = dist0_[bi][bj];
            d_vector = (ca_xyz_[bi]-ca_xyz_[bj]+dxyz[bi]-dxyz[bj]);
            d_new = d_vector.length();
            correction = d_vector*( (d_new-d0)/d_new/2.0 );
            correct_xyz_[bi] -= correction*sf_;
            correct_xyz_[bj] += correction*sf_;
            cws_[bi] += sf_;
            cws_[bj] += sf_;
          }
        }

        int np = indices.size();
        for(n=0;n<np;n++)
        {
          bi = indices[n][0];
          bj = indices[n][1];
          d0 = dist0_[bi][bj];
          if( d0 == 0 ) continue;
          d_vector = (ca_xyz_[bi]-ca_xyz_[bj]+dxyz[bi]-dxyz[bj]);
          d_new = d_vector.length();
          correction = d_vector*( (d_new-d0)/d_new/2.0 );
          correct_xyz_[bi] += correction;
          correct_xyz_[bj] -= correction;
          cws_[bi] +=1.0;
          cws_[bj] +=1.0;
        }
        // add correction to dxyz //
        for(n=0;n<n_ca_;n++)
        {
          this_dxyz_[n] = dxyz[n] + ( correct_xyz_[n]/cws_[n] ) ;
        }

        return this_dxyz_;
      }


        scitbx::af::shared <FloatType> Histogram(scitbx::af::const_ref <FloatType> array, FloatType dMax, int n_slots)
        {
           scitbx::af::shared <FloatType> hist;
           scitbx::histogram <FloatType,int> h(array, 0.0, dMax, n_slots);
           hist.push_back((h.slots()[0]*0.75+h.slots()[1]*0.25)/total_);
           for(int i=1;i<n_slots-1;i++)
             hist.push_back(((h.slots()[i-1]+h.slots()[i+1])/4.0 + h.slots()[i]/2.0)/total_);
           hist.push_back((h.slots()[n_slots-2]*0.25+h.slots()[n_slots-1]*0.75)/total_);
           return hist;
        }

        void updateDistArray(scitbx::af::const_ref <scitbx::vec3 < FloatType > > crd)
        {
          int count=0;
          total_=static_cast<FloatType> (n_ca_*(n_ca_-1)/2);
          scitbx::vec3 <FloatType> Vdiff;
          for(int i=0;i<n_ca_-1;i++)
          {  for(int j=i+1;j<n_ca_;j++)
            {
                Vdiff = crd[ca_indx_[i]]-crd[ca_indx_[j]];
                dist_array[count]=sqrt(Vdiff*Vdiff);
                count++;
            }
          }
        }

        void updateDistArrayAll(scitbx::af::const_ref <scitbx::vec3 < FloatType > > crd)
        {
          int count=0;
          total_=static_cast<FloatType> (natom_*(natom_-1)/2);
          if(dist_array.size() < natom_*(natom_-1)/2) {
            dist_array.reserve(natom_*(natom_-1)/2);
            for(int i=0;i<natom_*(natom_-1)/2-n_ca_*(n_ca_-1)/2;i++)
                dist_array.push_back(0);
                }
          scitbx::vec3 <FloatType> Vdiff;
          for(int i=0;i<natom_-1;i++)
          {  for(int j=i+1;j<natom_;j++)
            {
                Vdiff = crd[i]-crd[j];
                dist_array[count]=sqrt(Vdiff*Vdiff);
                count++;
            }
          }
        }

        scitbx::af::shared <FloatType> getDistArray()
        {
                return dist_array;
        }

        scitbx::af::shared<scitbx::vec3< FloatType> > ca_mode( int vec_indx)
        {
          scitbx::af::shared< scitbx::vec3< FloatType > > this_mode( n_ca_, scitbx::vec3<FloatType>(0.0) );
          FloatType upper(mean_disp_*4), lower((-1.0)*upper); //set boundary
          int indx(0), j;
          for(int i=0;i<n_ca_;i++)
          {
            for(j=0;j<3;j++)
            {
              this_mode[i][j] = eigenvec[indx][vec_indx];
              indx++;
              if( this_mode[i][j] > upper ) this_mode[i][j] = upper;
              if( this_mode[i][j] < lower ) this_mode[i][j] = lower;
            }
          }
          return this_mode;
        }


        scitbx::af::shared<scitbx::vec3< FloatType > >
        project2all( scitbx::af::const_ref< scitbx::vec3<FloatType> > ca_dxyz )
        {
          scitbx::af::shared <scitbx::vec3< FloatType > > Proj(natom_, scitbx::vec3<FloatType>(0.0) );
          for(int j=0;j<ca_indx_[1];j++)
            Proj[j] = ca_dxyz[1];
          for(int i=1;i<n_ca_-2;i++)
            for( int j=ca_indx_[i];j<ca_indx_[i+1];j++ )
              Proj[j] = ca_dxyz[i];
          for(int j=ca_indx_[n_ca_-2];j<natom_;j++)
            Proj[j] = ca_dxyz[n_ca_-2];
          return Proj;
        }

        scitbx::af::shared<scitbx::vec3< FloatType > >Project2All(int vec_indx)
        {
          scitbx::af::shared <scitbx::vec3< FloatType > > Proj(natom_, scitbx::vec3<FloatType>(0.0) );
          //int back_indx(eigenvec.dim1()-vec_indx);
          int back_indx(vec_indx);
          Damp_vector(back_indx);
          for(int j=0;j<ca_indx_[1];j++)
            for(int x_i=0;x_i<3;x_i++)
              Proj[j][x_i] = (eigenvec[3+x_i][back_indx]);
          for(int i=1;i<n_ca_-2;i++)
          {
            for(int j=ca_indx_[i];j<ca_indx_[i+1];j++)
              for(int x_i=0;x_i<3; x_i++)
                Proj[j][x_i] = (eigenvec[i*3+x_i][back_indx]);
          }
          int ca_last=ca_indx_[n_ca_-2];
          for(int j=ca_last;j<natom_;j++)
             for(int x_i=0; x_i<3; x_i++)
               Proj[j][x_i] = (eigenvec[(n_ca_-2)*3+x_i][back_indx]);
          return Proj;
        }

        void Damp_vector(int indx)
        {
         FloatType two_mean(mean_disp_*2), n_two_mean((-1.0)*two_mean);

         int dimen = eigenvec.dim1();
         for(int i=0;i<dimen;i++)
         {
           if(eigenvec[i][indx] > two_mean)
           {
                eigenvec[i][indx] = two_mean;
           }
           else if(eigenvec[i][indx] < n_two_mean)
           {
                eigenvec[i][indx] = n_two_mean;
           }
         }
        }

        private:
          FloatType cutoff_, mean_disp_, sf_;
          int natom_, n_ca_;
          FloatType total_;
          scitbx::af::shared <int> ca_indx_;
          scitbx::af::shared <scitbx::vec3<FloatType> > ca_xyz_;
          scitbx::af::shared <scitbx::vec3<FloatType> > this_dxyz_;
          scitbx::af::shared <scitbx::vec3<FloatType> > correct_xyz_;
          scitbx::af::shared <FloatType > cws_;
          TNT::Array2D <FloatType> Hessian;
          TNT::Array2D <FloatType> dist0_;
          TNT::Array1D <FloatType> eigenval;
          scitbx::af::shared <FloatType> dist_array;
          TNT::Array2D <FloatType> eigenvec;
   }; //end of elastic class


} //namepsace refine
} //namespace sastbx
