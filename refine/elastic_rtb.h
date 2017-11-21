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
#include <tntbx/include/jama_qr.h>
#include <vector>

namespace sastbx { namespace refine {

   template<typename FloatType>
   class ProjMat
   {
    //Sparse matrix for RTB projection
     public:
       ProjMat() {nrows_ = ncols_ = 0; }
       ~ProjMat() {}
     private:
       std::vector< TNT::Array2D<FloatType> > D_;
       scitbx::af::shared<int> start_na_;
       int nrows_,ncols_;
     public:
       void push_back(TNT::Array2D <FloatType> P)
       {
         D_.push_back(P);
         start_na_.push_back(nrows_);
         nrows_ += P.dim1();
         ncols_ += P.dim2();
       }
       int  num_rows(){return nrows_;}
       int  num_cols(){return ncols_;}

       FloatType operator()(int const& i, int const& j)
       {
         SCITBX_ASSERT(i <= nrows_);
         SCITBX_ASSERT(j <= ncols_);
         int nb(j/6); //block number from j
         int j_in_block(j-nb*6); //offset
         int i_in_block(i-start_na_[nb]);
         if (i_in_block < 0) return 0;
         if (i_in_block < D_[nb].dim1() &&
             j_in_block < D_[nb].dim2())
           return D_[nb][i_in_block][j_in_block];
         return 0;
       }
       void clear() { start_na_.clear(); D_.clear(); }
    }; //end of ProjMat


   template<typename FloatType>
   class elastic_rtb
   {
        public:
          elastic_rtb(scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& xyz,
                  scitbx::af::const_ref <int> const& block_start,
                  FloatType const& cutoff,
                  FloatType const& sf):
          cutoff_(cutoff), sf_(sf), cutoff2_(cutoff*cutoff),
          MAX_N_BLOCK_(1000), MIN_N_ATOM_(10)
          {
           natom_ = xyz.size();
           n_block_  = block_start.size();

           for(int i=0;i<n_block_;i++)
             block_start_.push_back(block_start[i]);
           make_blocks();
           build_proj_mat(xyz); //Projection
           Hessian = calc_hessian(xyz); //all-atom hessian
           HessianP = calc_hessian_p(Hessian); //projected PHP
           JAMA::Eigenvalue <FloatType> eig(HessianP);
           eig.getV(eigenvec);
           eig.getRealEigenvalues(eigenval);
           //std::cout<<"E-values"<<std::endl;
           //std::cout << eigenval;
           //std::cout<<"E-vectors"<<std::endl;
           //for(int i=0;i<n_block_*3;i++)
           //std::cout <<eigenvec[i][0]<<" "<<eigenvec[i][7]<<std::endl;
//           eig_vec_all = deproject();

        }

        int nmode(FloatType percentage, int n_max_mode)
        {
          FloatType total(0.0), partial(0.0);
          int n_val( eigenval.dim() );
          int ii;
          for(ii=0;ii<n_val;ii++)
            total += eigenval[ii];
          total *= percentage;
          for(ii=0;ii<n_max_mode;ii++)
          {
            partial += eigenval[ii];
            if( partial >= total) break;
          }
          return ii;
        }

        void make_blocks()
        {
/*          if(natom_/MIN_N_ATOM_ > MAX_N_BLOCK_) natom_in_block=natom_/MAX_N_BLOCK_;
          else natom_in_block=MIN_N_ATOM_;

          n_block_ = int(natom_/natom_in_block);
*/
          for(int i=0;i<n_block_-1;i++)
          {
            natom_in_block_.push_back(block_start_[i+1]-block_start_[i]);
//            block_start_.push_back(natom_in_block*i);
          }
          natom_in_block_.push_back(natom_-block_start_[n_block_-1]);
          // if there is any extra atoms left, add them to the last block
//          natom_in_block_[n_block_-1] += (natom_-natom_in_block*n_block);
        }

        void build_proj_mat(scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& xyz)
        {
          for (int nb = 0; nb < n_block_; nb++)
          {
            TNT::Array2D<FloatType> P(3*natom_in_block_[nb],6, 0.0);
            scitbx::vec3<FloatType> CoM(0,0,0);
            for (int nj = 0; nj < natom_in_block_[nb]; nj++)
              CoM += xyz[block_start_[nb]+nj];
            CoM /= static_cast<FloatType>(natom_in_block_[nb]);
            for (int nj = 0; nj < natom_in_block_[nb]; nj++)
            {
              scitbx::vec3<FloatType> rdiff(xyz[block_start_[nb]+nj]-CoM);
              P[3*nj  ][0] = 1.0;
              P[3*nj+1][1] = 1.0;
              P[3*nj+2][2] = 1.0;
              P[3*nj  ][4] = rdiff[2];
              P[3*nj  ][5] =-rdiff[1];
              P[3*nj+1][5] = rdiff[0];
              P[3*nj+1][3] = -P[3*nj][4];
              P[3*nj+2][3] = -P[3*nj][5];
              P[3*nj+2][4] = -P[3*nj+1][5];
            }
            JAMA::QR<FloatType> PQR(P);
            ProjMat_.push_back(PQR.getQ());
          }
        }

        TNT::Array2D<FloatType> calc_hessian(scitbx::af::const_ref <scitbx::vec3 <FloatType> > const& xyz)
        {
           Hessian=TNT::Array2D<FloatType>(natom_*3,natom_*3,0.0);

           int hi,hj,indx_i,indx_j;
           scitbx::vec3 <FloatType> vecDiff;
           FloatType distsqr, Hessij;
           FloatType scale_factor;

           for(int i=0;i<natom_;i++)
           { hi = i*3;
             for(int j=i+1;j<natom_;j++)
             {  hj=j*3;
                vecDiff = xyz[i]-xyz[j];
                distsqr = vecDiff * vecDiff;
                if( distsqr == 0 )
                { std::cout<<"distance between atoms is zero!("
                           <<i<<","<<j<<")"<<std::endl; }
                if(j-i < 4) { scale_factor = sf_;}
                else { scale_factor = 1;}
                if(distsqr <= cutoff2_)
                {
                  for(int ii=0;ii<3;ii++)
                    for(int jj=ii;jj<3;jj++)
                    {
                        Hessij = vecDiff[ii]*vecDiff[jj]/distsqr*scale_factor;
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
           }//end of for-loop natom_i_
           return Hessian;
        }

        TNT::Array2D<FloatType> calc_hessian_p(TNT::Array2D<FloatType> HA) //Hp = P*H*P
        {
          SCITBX_ASSERT(HA.dim2() == ProjMat_.num_rows());
          int I = HA.dim1();
          int J = ProjMat_.num_cols();
          TNT::Array2D<FloatType> HAP(I,J,0.0);
          for (int i=0; i<I; i++)
          {
            for (int j=0; j<J; j++)
            {
              int nb(j/6);
              int mink(block_start_[nb]*3);
              int maxk(block_start_[nb]*3+natom_in_block_[nb]*3);
//              std::cout<<nb<<" "<<mink<<" "<<maxk<<std::endl;
              for (int k=mink; k<maxk; k++)
                HAP[i][j] += HA[i][k] * ProjMat_(k,j);
            }
          }

          SCITBX_ASSERT(ProjMat_.num_rows() == HAP.dim1() ); //dim1 == # of rows
          I = ProjMat_.num_cols(); //transpose
          J = HAP.dim2();
          //std::cout<<I<<" "<<J<<" hessianp"<<std::endl;
          TNT::Array2D<FloatType> HessianP(I,J,0.0);

          for (int i=0; i<I; i++)
          {
            for (int j=i; j<J; j++)
            {
              int nb(i/6);
              int mink(block_start_[nb]*3);
              int maxk(block_start_[nb]*3+natom_in_block_[nb]*3);
              for (int k=mink; k<maxk; k++)
                HessianP[i][j] += ProjMat_(k,i) * HAP[k][j];
              HessianP[j][i] = HessianP[i][j];
            }
          }
          return HessianP;
        }


        TNT::Array2D<FloatType> deproject() //Vec=P*eigvec
        {
          SCITBX_ASSERT(eigenvec.dim1() == ProjMat_.num_cols());
          int I = ProjMat_.num_rows();
          int J = eigenvec.dim2();
          TNT::Array2D<FloatType> eig_all(I,J,0.0);
          int nb(0);
          for (int i=0; i<I; i++)
          {
            int ii(block_start_[nb]*3+natom_in_block_[nb]*3);
            if (i > ii) nb++;
            int mink(6*nb);
            int maxk(6*(nb+1));
            for (int j=0; j<J; j++)
            {
              for (int k=mink; k<maxk; k++)
                eig_all[i][j] += ProjMat_(i,k) * eigenvec[k][j];
            }
          }
          return eig_all;
        }

        scitbx::af::shared< scitbx::vec3< FloatType > > mode(int n)
        {
          SCITBX_ASSERT(eigenvec.dim1() == ProjMat_.num_cols());
          int I = ProjMat_.num_rows(), nb(0);
          scitbx::af::shared<FloatType> mode_n(I,0.0);
          scitbx::af::shared< scitbx::vec3< FloatType > > this_mode(natom_,scitbx::vec3<FloatType>(0.0) );
          for(int i=0;i<I;i++)
          {
            int ii(block_start_[nb]*3+natom_in_block_[nb]*3);
            if (i > ii) nb++;
            int mink(6*nb);
            int maxk(6*(nb+1));
            for (int k=mink; k<maxk; k++)
              mode_n[i] += ProjMat_(i,k) * eigenvec[k][n];
          }
          int indx(0), j;
          for(int i=0; i<natom_; i++)
            for(j=0;j<3;j++)
            {
              this_mode[i][j] = mode_n[indx];
              indx++;
            }
          return this_mode;
        }

        scitbx::af::shared<scitbx::vec3< FloatType > >
        project2all( scitbx::af::const_ref< scitbx::vec3<FloatType> > ca_dxyz )
        {
          int I = ProjMat_.num_rows(), nb(0);
          scitbx::af::shared<FloatType> mode_n(I,0.0);
          scitbx::af::const_ref<FloatType> dx(ca_dxyz);
          for(int i=0;i<I;i++)
          {
            int ii(block_start_[nb]*3+natom_in_block_[nb]*3);
            if (i > ii) nb++;
            int mink(6*nb);
            int maxk(6*(nb+1));
            for (int k=mink; k<maxk; k++)
              mode_n[i] += ProjMat_(i,k) * dx[k];
          }
          scitbx::af::shared <scitbx::vec3< FloatType > > Proj(natom_, scitbx::vec3<FloatType>(mode_n) );
/*
          int indx(0), j;
          for(int i=0; i<natom_; i++)
            for(j=0;j<3;j++)
            {
              Proj[i][j] = mode_n[indx];
              indx++;
            }
*/
          return Proj;
        }

        private:
          FloatType cutoff_, mean_disp_, sf_, cutoff2_;
          int natom_, n_block_, MAX_N_BLOCK_, MIN_N_ATOM_;
          FloatType total_;
          scitbx::af::shared <int> block_start_;
          scitbx::af::shared <int> natom_in_block_;
          TNT::Array2D <FloatType> Hessian;
          TNT::Array2D <FloatType> HessianP;
          TNT::Array1D <FloatType> eigenval;
          TNT::Array2D <FloatType> eigenvec;
          TNT::Array2D <FloatType> eig_vec_all;
          ProjMat<FloatType> ProjMat_;
   }; //end of elastic RTB class




} //namepsace refine
} //namespace sastbx
