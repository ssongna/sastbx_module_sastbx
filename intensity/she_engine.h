//! Peter Zwart April 05, 2005
#ifndef SASTBX_INTENSITY_SHE_H
#define SASTBX_INTENSITY_SHE_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <scitbx/array_family/shared.h>
#include <map>

#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <scitbx/wigner/wigner3j.h>
#include <string>
#include <iomanip>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <complex>

#include <sastbx/intensity/integrated.h>
#include <sastbx/intensity/besselR.h>
#include <sastbx/intensity/model.h>
#include <time.h>
//#include <cctbx/math/cos_sin_table.h>

using boost::math::spherical_harmonic; //special_functions::spherical_harmonic;
using boost::math::sph_bessel;
using boost::math::legendre_p;
using scitbx::wigner::wigner3j_fast;
using sastbx::intensity::model;
//using sastbx::rigid_body::wigner3j_fast;
using namespace std;

namespace sastbx { namespace intensity {

  template <typename FloatType>
  class she_engine
  {
    public:
      she_engine(
                     model <FloatType> const& model,
                     scattering_library< FloatType > scat_lib,
                     int const& max_i,
                     int const& max_L,
                     FloatType const& f_resl,
                     FloatType const& q_step,
                     FloatType const& max_z,
                     FloatType const& delta,
                     FloatType const& rho,
                     FloatType const& d_rho
                  ):
      model_( model ),
      scat_lib_( scat_lib),
      center_(0,0,0),
      max_i_(max_i), max_L_(max_L), delta_(delta),rho_(rho), d_rho_(d_rho), f_w_max_(0.0), f_w_min_(10000.0),
      f_resl_(f_resl),q_step_(q_step), max_z_(max_z), wigner_fast_(max_L_*3+1), add_wigner(true)  //bessel_array_(0.1,max_z,max_L)
      {
        n_=model_.size();
        for(int i=0; i < n_; i++)
        {
         xyz_.push_back(model_.get_xyz(i));
         rtp_.push_back(model_.get_rtp(i));
         radius_.push_back(model_.get_radius(i));
         atom_cs_t_.push_back(scitbx::vec2<FloatType>(cos(rtp_[i][1]), sin(rtp_[i][1])) );
        }
        Grid();
        F_omega_array();
        Build_integration_map();
        int n_q = scat_lib_.n_q();
        I_all_.reserve(n_q);
        I_A_.reserve(n_q);
        I_B_.reserve(n_q);
        I_C_.reserve(n_q);
        Compute_Coefs();
      }//end of she_engine constructor

// calculate I(s) using bessel spherical Harmonic expansion formula
      void Compute_Coefs()
      {
        int n_q = scat_lib_.n_q();
        FloatType q, this_I, tmp_float;
        scitbx::af::shared< std::complex<FloatType> > sph_array;
        scitbx::af::shared< std::complex<FloatType> > sph_B;
        scitbx::af::shared< std::complex<FloatType> > i_pow;
        scitbx::af::shared< FloatType > J_l_array;
        scitbx::af::shared< FloatType > int_l;
        std::complex<FloatType> complex_I(0,1), co_l, tmp_complex,tmp_dummy_c, zero(0,0), i_pow_l;
        std::complex<FloatType> sph,Integral_r, sph_b;
        int accum, sph_count, lm_size;
        FloatType four_pi=scitbx::constants::pi * 4.0;
        int max_grid = fibonacci(max_i_);
        FloatType integration_unit=four_pi/static_cast<FloatType>(max_grid+1);
        FloatType ia,ib,ic;
        //time_t t;

        A_array_.clear();
        B_array_.clear();
        C_array_.clear();
        lm_size = (max_L_+2)*(max_L_+1)/2;
        J_l_array.reserve(n_);
        int_l.reserve(max_grid+1); //bessel integral at SPH order l (max_grid+1) points

        //caculate A_lm and C_lm
        for(int l=0;l<=max_L_;l++)
        {
         i_pow.push_back( pow(complex_I, l) );
         for(int m=0;m<=l;m++) // Y_l^m* = (-1)^m Y_l^(-m)
         {
          for(int i=0;i<n_;i++)
          {
                sph_array.push_back(conj(spherical_harmonic(l,m,rtp_[i][1],rtp_[i][2])));
          }
          for(int i=0;i<=max_grid;i++)
          {
                sph_B.push_back(conj(spherical_harmonic(l,m,grid_[i][0],grid_[i][1])));
          }
         }
        } //end of l

        for (int index=0;index<n_q;index++)
        {
        scitbx::af::shared< std::complex<FloatType> > A_q;
        scitbx::af::shared< std::complex<FloatType> > B_q;
        scitbx::af::shared< std::complex<FloatType> > C_q;
        A_q.reserve(lm_size);
        B_q.reserve(lm_size);
        C_q.reserve(lm_size);

         accum=0;
         sph_count = 0;
         q  = scat_lib_.q(index);
        for(int l=0;l<=max_L_; l++)
        {
         for(int j=0;j<n_;j++) J_l_array[j]=sph_bessel(l,q*rtp_[j][0]); //bessel_array_.get_value(l,q*rtp_[j][0]);
         i_pow_l = i_pow[l];
         for(int m=0;m<=l;m++)
         {
          tmp_complex=zero;
          tmp_dummy_c=zero;
          for(int j=0;j<n_;j++)
          {
                sph_b=J_l_array[j]*sph_array[sph_count++];
                tmp_complex += model_.sf_array_[index][j]*sph_b;
                tmp_dummy_c += model_.dummy_sf_array_[index][j]*sph_b;
          } //end of j

          A_q[accum]= tmp_complex*i_pow_l;
          C_q[accum]= tmp_dummy_c*i_pow_l;  //increment index of A,C arrays, ONCE
          accum++;
         }//end of m
        }//end of l

//      Calculate B_Array
        accum=0;
        sph_count=0;
        for(int l=0;l<=max_L_;l++)
        {
         i_pow_l = i_pow[l];
         for(int ii=0;ii<=max_grid;ii++)
           int_l[ii]=Integral_[l].Get_integral(F_array_[ii], q);

         for(int m=0;m<=l;m++)
         {
          Integral_r = 0.0;
          for(int ii=0;ii<=max_grid;ii++)
          {
                Integral_r += (sph_B[sph_count]*int_l[ii]);
                sph_count++;
          }
          B_q[accum]=Integral_r*integration_unit*i_pow_l;
          accum++;
         } //end of m for B_q_array
        }//end of l for B__q_array

        A_array_.push_back(A_q);
        B_array_.push_back(B_q);
        C_array_.push_back(C_q);
      }//end of q
    }

     scitbx::af::shared< FloatType > I()
     {
        int n_q = scat_lib_.n_q();
        scitbx::af::shared <FloatType> scales;
        for (int index=0; index < n_q; index++)
           scales.push_back(1.0);
        return Iscale(scales.const_ref());
     }

     scitbx::af::shared< FloatType > Iscale(scitbx::af::const_ref<FloatType> scales)
     {
        int n_q = scat_lib_.n_q();
        FloatType q, scale,this_I,ia,ib,ic;
        int accum;
        FloatType four_pi=scitbx::constants::pi * 4.0;
        I_all_.clear();
        I_A_.clear();
        I_B_.clear();
        I_C_.clear();

        for (int index=0;index<n_q;index++)
        {
         q  = scat_lib_.q(index);
         scale = scales[index];

        this_I=0.0;
        ia=0.0;
        ib=0.0;
        ic=0.0;

        accum=0;
        for(int l=0;l<=max_L_;l++)
        { // Al0 terms should be scaled to half, then the whole thing will be scaled back
          ia += norm(A_array_[index][accum])/2.0;
          ic += norm(C_array_[index][accum])/2.0;
          ib += norm(B_array_[index][accum])/2.0;
          this_I += norm(A_array_[index][accum]-rho_*scale*C_array_[index][accum]+d_rho_*B_array_[index][accum])/2.0;
          accum++;

         for(int m=1;m<=l;m++)
         {
                ia += norm(A_array_[index][accum]);
                ic += norm(C_array_[index][accum]);
                ib += norm(B_array_[index][accum]);
                this_I += norm(A_array_[index][accum]-rho_*scale*C_array_[index][accum]+d_rho_*B_array_[index][accum]);
                accum++;
         }
        }
         //cout<<setw(15)<<q<<setw(15)<<four_pi*this_I<<setw(15)<<four_pi*ia<<setw(15)<<four_pi*ib<<setw(15)<<four_pi*ic<<endl;
         I_all_.push_back( four_pi*this_I *2.0);
         I_A_.push_back(four_pi*ia *2.0);
         I_C_.push_back(four_pi*ic *2.0*pow(rho_*scale,2.0));
         I_B_.push_back(four_pi*ib *2.0*pow(d_rho_,2.0));
        }//end of I(q)
        return (I_all_);
      }

      scitbx::af::shared < FloatType > get_IA() {return I_A_;}
      scitbx::af::shared < FloatType > get_IB() {return I_B_;}
      scitbx::af::shared < FloatType > get_IC() {return I_C_;}

      void UpdateSolventParams(FloatType rho, FloatType drho)
      {
        rho_=rho;
        d_rho_=drho;
      }

     // Build the quasiuniform grid
      void Grid()
      {
        int i,max_grid,second_max;
        FloatType theta, phi, f_max, two_pi;

        max_grid=fibonacci(max_i_);
        second_max=fibonacci(max_i_-1);
        f_max = static_cast<FloatType> (max_grid);
        two_pi = scitbx::constants::pi * 2.0;

        for(i=0;i<=max_grid;i++)
         {
                scitbx::vec2<FloatType> omega;
                theta =acos(1.0-2.0*i/f_max); //[0,pi]
                phi = two_pi * ((i*max_grid) % second_max)/static_cast<FloatType>(second_max);
                omega[0]=theta;
                omega[1]=phi;
                grid_.push_back(omega);
                grid_cs_t_.push_back(scitbx::vec2<FloatType>( cos(theta), sin(theta) ) );
         }
        return ;
      }//end of Griding

      void F_omega_array()
      {
        int i, total, max_grid, indx, second_max,span;

        total=grid_.size();
        F_array_.clear();
        F_array_.reserve(total);
        for(i=0;i<total;i++) F_array_[i]=0.0;

        max_grid=total-1;
        second_max=fibonacci(max_i_-1);

        for(i=0;i<n_;i++)
        {
          for(int j=0;j<total;j++)
          {
                update_F_w(j,i);
          }
        }
        for(int j=1;j<total-1;j++)
        { if(F_array_[j]==0.0) F_array_[j]=(F_array_[j-1]+F_array_[j+1])*0.5;}
        if(F_array_[0]==0.0) F_array_[0]=F_array_[1];
        if(F_array_[total-1]==0.0) F_array_[total-1]=F_array_[total-2];

        for(int j=0;j<total;j++) {
          if(f_w_min_> F_array_[j]) f_w_min_=F_array_[j];
          if(f_w_max_< F_array_[j]) f_w_max_=F_array_[j];
          //std::cout<<F_array_[j]<<std::endl;
        }
        return ;
      }

      void update_F_w(int direction_j, int atom_i)
      {
        ct0=grid_cs_t_[direction_j][0]; //cosine(t0)
        st0=grid_cs_t_[direction_j][1]; //sine(t0)
        p0=grid_[direction_j][1];

        r1=rtp_[atom_i][0];
        t1=rtp_[atom_i][1];
        p1=rtp_[atom_i][2];
        cos_ad = ct0*atom_cs_t_[atom_i][0]+st0*atom_cs_t_[atom_i][1]*cos(p1-p0);
        if(cos_ad >=0 || n_==1)
        {
          sin_ad = sqrt(1.0 - cos_ad*cos_ad);
          distance = r1*sin_ad;
          if(distance < (1.5+radius_[atom_i])) //need to be corrected to r_w + r_[i]
          {
             distance = r1*cos_ad;
             if(F_array_[direction_j] < distance+radius_[atom_i]*0.5)
             { F_array_[direction_j] = distance + radius_[atom_i]*0.5;}
          }
        }
        return;
      }
      scitbx::vec3<FloatType> Area_Volume()
      {
        FloatType AREA=0.0, Volume_IN=0.0, Volume_OUT=0.0, PI=scitbx::constants::pi, two_PI;
        int max_grid=fibonacci(max_i_);
        FloatType np = static_cast<FloatType> (max_grid+1);
        two_PI= 2.0*PI;
        for(int i=0;i<=max_grid;i++)
        {
         AREA += pow(F_array_[i]-0.3,2.0);
         Volume_IN += pow(F_array_[i],3.0);
         Volume_OUT += pow(F_array_[i]+delta_,3.0);
        }
        AREA *= (PI*4.0/np);
        Volume_IN *= (PI*4.0/np/3.0);
        Volume_OUT *= (PI*4.0/np/3.0);
        scitbx::vec3<FloatType> area_volume_shell(AREA, Volume_IN, Volume_OUT-Volume_IN);
        return area_volume_shell;
      }

      scitbx::vec3<FloatType> Area_Volume2()
      {
        FloatType AREA=0.0, Volume_IN=0.0, Volume_OUT=0.0, PI=scitbx::constants::pi, two_PI;
        int max_grid=fibonacci(max_i_);
        FloatType np = static_cast<FloatType> (max_grid+1), unit;
        two_PI= 2.0*PI;

        Find_neighbors();
        for(int i=0;i<max_grid-2;i++)
        {
         unit=quadra(i);
         AREA += unit;
         Volume_IN += pow(F_array_[i],3.0)*unit;
         Volume_OUT += pow(F_array_[i]+delta_,3.0)*unit;
        }
        Volume_IN /= (3.0);
        Volume_OUT /= (3.0);
        scitbx::vec3<FloatType> area_volume_shell(AREA, Volume_IN, Volume_OUT-Volume_IN);
        return area_volume_shell;
      }

      void Find_neighbors()
      {
       FloatType distij,min1,min2,min3;
       int m1_indx, m2_indx, m3_indx;
       int max_grid=fibonacci(max_i_);

       for(int i=0;i<max_grid-2;i++)
        {
          min1=min2=min3=100000000;
          m1_indx=m2_indx=m3_indx=0;
          for(int j=i+1;j<=max_grid;j++)
          {
                distij=(grid_[i]-grid_[j]).length_sq();
                if(distij < min1)
                {
                  min3=min2;
                  min2=min1;
                  min1=distij;
                  m3_indx=m2_indx;
                  m2_indx=m1_indx;
                  m1_indx=j;
                }
                else if(distij < min2)
                {
                  min3=min2;
                  min2=distij;
                  m3_indx=m2_indx;
                  m2_indx=j;
                }
                else if(distij < min3)
                {
                  min3=distij;
                  m3_indx=j;
                }
          }
        scitbx::vec3<int> neighbor_i(m1_indx, m2_indx, m3_indx);
        neighbors_.push_back(neighbor_i);
        }
        return;
      }

      FloatType quadra(int i)
      {
        int j;
        FloatType x,y,z,theta,phi;
        theta=grid_[i][0]; phi=grid_[i][1];
        x=sin(theta);
        y=x*sin(phi);
        x=x*cos(phi);
        z=cos(theta);
        scitbx::vec3<FloatType> V0(x,y,z);
        V0 *= F_array_[i];
        j = neighbors_[i][1];
        theta=grid_[j][0]; phi=grid_[j][1];
        x=sin(theta);
        y=x*sin(phi);
        x=x*cos(phi);
        z=cos(theta);
        scitbx::vec3<FloatType> V1(x,y,z);
        V1 *= F_array_[j];
        j = neighbors_[i][0];
        theta=grid_[j][0]; phi=grid_[j][1];
        x=sin(theta);
        y=x*sin(phi);
        x=x*cos(phi);
        z=cos(theta);
        scitbx::vec3<FloatType> V2(x,y,z);
        V2 *= F_array_[j];
        j = neighbors_[i][2];
        theta=grid_[j][0]; phi=grid_[j][1];
        x=sin(theta);
        y=x*sin(phi);
        x=x*cos(phi);
        z=cos(theta);
        scitbx::vec3<FloatType> V3(x,y,z);
        V3 *= F_array_[j];

        return (fabs((V0-V1)*(V1-V2))+abs((V1-V2)*(V2-V3))+abs((V2-V3)*(V3-V0))+abs((V3-V0)*(V0-V1)))/8.0;
      }

      int fibonacci(int n)
      {
        if(n==0) {return 0;}
        else if(n == 1) { return 1;}
        else { return fibonacci(n-1) + fibonacci(n-2);}
      }

      void Build_integration_map()
      {
        FloatType q_max;
        q_max=scat_lib_.q(scat_lib_.n_q()-1);
        upper_limit_= f_w_max_+f_resl_;
        lower_limit_=std::max(0.0,f_w_min_-f_resl_);
        for(int l=0;l<=max_L_;l++)
        {
                integrated <FloatType> it(lower_limit_,upper_limit_,delta_,f_resl_,q_step_,q_max,l);
                Integral_.push_back(it);
        }
        return;
      }

      void update_coordinate(scitbx::af::const_ref <scitbx::vec3 <FloatType> > xyz, scitbx::af::const_ref <int> indx)
      {
        //print_coordinate(indx);
        int atom_i, n_atom = indx.size();
        int max_grid = fibonacci(max_i_);
        scitbx::vec3 <FloatType> this_center(0,0,0);
        scitbx::vec3 <FloatType> old_center(0,0,0);

        FloatType this_max=0, this_min=10000;
        if(n_atom < n_/2)
        {
         for(int i=0;i<n_atom;i++)
          { atom_i=indx[i];
          old_center += xyz_[atom_i];
          this_center += xyz[i];
          }
          old_center = old_center/n_atom;
          this_center = this_center/n_atom;
          center_ = this_center-old_center;
        }
        else
        {
          for(int i=0;i<n_atom;i++)
          {
            xyz_[indx[i]]=xyz[i];
          }
          for(int i=0;i<n_;i++) this_center += xyz_[i];
        center_ = this_center/n_;
        }
        for(int i=0;i<n_;i++)
        { xyz_[i] -= center_;
          rtp_[i]=polar(xyz_[i], i);
          if(rtp_[i][0] > this_max) this_max=rtp_[i][0];
        }
/*        if((this_max*(scat_lib_.q(scat_lib_.n_q()-1)+0.01)) > max_z_)
        {       max_z_ = this_max*(scat_lib_.q(scat_lib_.n_q()-1)+0.02);
                besselr<FloatType> b(0.1,max_z_,max_L_);
                bessel_array_=b;
                std::cout<<"Bessel_array is updated with max_z = "<<max_z_<<std::endl;
        }
 */
       F_omega_array(); //updated f_w_max_ and f_w_min_ as well

        if(f_w_max_ > upper_limit_ || f_w_min_<lower_limit_)
        {
          Integral_.clear();
          Build_integration_map();
          std::cout<<"Integration map is updated with boundary = "<<lower_limit_<<":"<<upper_limit_<<std::endl;
        }

        Compute_Coefs();
        add_wigner = true;
        coef_array_.clear();
        //print_coordinate(indx);
        return;
      }

      scitbx::af::shared< FloatType >  get_spatial_correlation()
      { scitbx::af::shared< FloatType > correlation;
        int n1,n2,i,j,k;
        n1 = spatial_correlation_.size();
        n2 = scat_lib_.n_q();
        for(i=0; i<n1; i++)
          for(j=0; j<n2; j++)
            for(k=0; k<=j; k++)
             correlation.push_back( spatial_correlation_[i][j][k] );
        return correlation;
      }

      scitbx::af::shared<FloatType> get_expansion_coef(int q_index)
      {
        return expansion_coef_array_[q_index];
      }

      scitbx::af::shared<FloatType> get_all_coefs()
      {
        int n_q = scat_lib_.n_q();
        scitbx::af::shared<FloatType> B;
        for(int i=0;i<n_q;i++)
          B.insert( B.end(), expansion_coef_array_[i].begin(),expansion_coef_array_[i].end() );
        return B;
      }

      scitbx::af::shared<FloatType> expansion_coef(scitbx::af::const_ref< std::complex<FloatType> > I_lm_q )
      {
        scitbx::af::shared< FloatType > B_q;
        FloatType B_q_l;
        int count = 0;
        for(int l=0;l<=max_L_;l+=2)
        {
          B_q_l = 0.0;
          for(int m=-l; m<=l; m++)
          {
            B_q_l += norm( I_lm_q[count] );
            count++;
          }
          B_q.push_back(B_q_l);
        }
        return B_q;
      }



      void calc_spatial_correlation( scitbx::af::const_ref<FloatType> phi_array)
      {
        int n_total = (max_L_+1)*(max_L_+1);
        int n_angle = phi_array.size();
        int n_q = scat_lib_.n_q();
        int i, index, counter;
        S_array_.clear();
        I_lm_array_.clear();
        time_t t;

        for(int index=0;index<n_q; index++)
        {
          scitbx::af::shared< std::complex<FloatType> > S_q;
          scitbx::af::shared< std::complex<FloatType> > I_lm_q;
        //  for(int i=0; i< n_total; i++)
        //    {
        //      S_q.push_back( A_array_[index][i] - rho_*C_array_[index][i] + d_rho_*B_array_[index][i]) ;
        //        S_q.push_back( A_array_[index][i] );
        //    }
          counter = 0;
          for(int l=0;l<=max_L_;l++)
          {
            counter += l; //fast forward to A_l_l; mirroring A_l^(-m)
            for(int m=l;m>0;m--)
            { 
              S_q.push_back( conj(A_array_[index][counter])*pow(-1.0,m) );
              counter--;
            }
            for(int m=0;m<=l;m++)
            { 
              S_q.push_back( A_array_[index][counter] );
              counter++;
            }
            
          }

          accum_ = 0;
          int count=0;
          for( int l=0; l<=max_L_; l+=2 )
          {
            for(int m=-l; m<=l; m++)
            {
              I_lm_q.push_back( calc_I_lm(l,m,S_q.const_ref() ) );
//            std::cout<<I_lm_q[count]<<" ";
              count++;
            }
//          std::cout<<"L= "<<l<<std::endl;
          }
    /***** calc expansion coefficients B_l *****/
          expansion_coef_array_.push_back( expansion_coef( I_lm_q.const_ref() ) );

          I_lm_array_.push_back(I_lm_q);
          add_wigner = false;
        }
    /******calc spatial correlation******/

        spatial_correlation_.clear();
        for( int angle=0; angle < n_angle; angle++)
        {
          spatial_correlation_.push_back( calc_correlation_phi( phi_array[angle] ) );
        }

        return;
      }

      std::complex< FloatType > calc_I_lm(int l, int m, scitbx::af::const_ref< std::complex< FloatType > > S_q)
      {
        std::complex< FloatType > I_lm_value(0,0), zero(0,0), sum_m1, sum_m2;
        FloatType coef_m1, coef_m2, wigner_l1, wigner_l2, wigner_lm1, wigner_lm2;
        FloatType four_pi=scitbx::constants::pi * 4.0,sqrt_value;
        int m2,l2, m1, l1;

        I_lm_value=zero;
// stupid way: loop around all the values of l1,l2,m1,m2
/*      for(l1=0;l1<=max_L_;l1++)
        {
          for(l2=0;l2<=max_L_;l2++)
          {
                sqrt_value=sqrt( (2*l1+1)*(2*l2+1)*(2*l+1)/four_pi );
                wigner_l1 = wigner_fast_.compute(l1,l2,l,0,0,0);
             for(m1=-l1;m1<=l1;m1++)
             {
                for(m2=-l2;m2<=l2;m2++)
                {
                  if(add_wigner)
                  {
                    coef_m1 = sqrt_value*wigner_l1*wigner_fast_.compute(l1,l2,l,m1,-m2,-m);
                    coef_array_.push_back(coef_m1);
                  }
                  else
                  {
                    coef_m1 = coef_array_[accum_++];
                  }
                  sum_m1 = coef_m1*S_q[l1*l1+l1+m1]*conj( S_q[l2*l2+l2+m2] );
                 if( m1 % 2 == 0)
                   I_lm_value += sum_m1;
                 else
                   I_lm_value -= sum_m1;

                }
             }
          }
        }
        if( m%2 == 0)
          return I_lm_value;
        else
          return -I_lm_value;


*/
//
        for( int l1=0;l1<=max_L_;l1++)
         for( int m1 = -l1; m1<=l1; m1++)
         {
          m2 = m1-m; //to satisfy m1-m2-m=0; i.e. m1+m2+m3=0 in standard repr
          for( int l2=abs(m2);l2<=max_L_;l2++)
          {
            if( (l1+l2+l) % 2 == 0)
            {
              if(add_wigner)
              {
                sqrt_value=sqrt( (2*l1+1)*(2*l2+1)*(2*l+1)/four_pi );
                wigner_l1 = wigner_fast_.compute(l1,l2,l,0,0,0);
                wigner_lm1 = wigner_fast_.compute(l1,l2,l,m1,-m2,-m);
                coef_m1 = wigner_l1* wigner_lm1* sqrt_value;
                coef_array_.push_back( coef_m1 );
                if( m != 0) {
                  wigner_l2 = wigner_l1; //wigner_fast_.compute(l2,l1,l,0,0,0);
                  wigner_lm2 = wigner_fast_.compute(l2,l1,l,m2,-m1,-m);
                  coef_m2 = wigner_l2* wigner_lm2* sqrt_value;
                  coef_array_.push_back( coef_m2 );
                  }
              }
              else
              {
                coef_m1 = coef_array_[accum_++];
                if( m != 0)
                  coef_m2 = coef_array_[accum_++];
              }
              sum_m1 = coef_m1*S_q[l1*l1+l1+m1]*conj( S_q[l2*l2+l2+m2] );
              if( m1 % 2 == 0)
                I_lm_value += sum_m1;
              else
                I_lm_value -= sum_m1;
              if( m != 0) {
                sum_m2 = coef_m2*S_q[l2*l2+l2+m2]*conj( S_q[l1*l1+l1+m1] );
                if( m2 % 2 == 0)
                  I_lm_value += sum_m2;
                else
                  I_lm_value -= sum_m2;
               }

             }
           }
          }
        if( m%2 == 0)
          return I_lm_value;
        else
          return -I_lm_value;

//      return I_lm_value*pow(-1,m);
      }


      scitbx::af::shared< scitbx::af::shared< FloatType > >  calc_correlation_phi( FloatType phi )
      {
        scitbx::af::shared< scitbx::af::shared< FloatType > > correlation;
        scitbx::af::shared< FloatType > p_l_cosine_phi;
        for(int i = 0; i <= max_L_; i ++)
        {
          p_l_cosine_phi.push_back( legendre_p( i, cos(phi) ) );
        }

        for( int i = 0 ; i < scat_lib_.n_q(); i++)
        {
           scitbx::af::shared< FloatType > corr_j;
           for( int j = 0 ; j <= i; j++ )
          {
                //corr_j.push_back( calc_corr_phi_i_j( p_l_cosine_phi.const_ref(), i, j ) );
                corr_j.push_back( calc_corr_phi_i_j( phi,i, j ) );
          }
          correlation.push_back( corr_j );
        }
        return correlation;
      }

//      FloatType calc_corr_phi_i_j( scitbx::af::const_ref<FloatType> p_l_cos_phi, int i, int j)
      FloatType calc_corr_phi_i_j( FloatType phi, int i, int j)
      {
        int accum = 0; // and the l starts with 1, because the zero term is ignored
        std::complex<FloatType> sum(0,0);
        std::complex<FloatType> sum_00(0,0);
        FloatType corr, cos_t_i2, sin_t_j2;

        for( int l = 0 ; l <= max_L_; l+=2 )
        {
          std::complex<FloatType> sum_l(0,0);
          for(int m = -l ; m <=l; m++)
          {
            sum_l += I_lm_array_[i][accum] * conj( I_lm_array_[j][accum] );
            accum++;
          }
      //    sum += p_l_cos_phi[l] *sum_l;
            sum_00 += sum_l;
            cos_t_i2 = pow( (scat_lib_.q(i)/12.56),2.0);
            sin_t_j2 = 1.0-cos_t_i2; //for same q, now
            sum += legendre_p(l, cos_t_i2+sin_t_j2*cos(phi))*sum_l;
        }
        //corr = abs( I_lm_array_[i][0]*conj(I_lm_array_[j][0]) );
        //comment out the following line to remove normalization
        //corr = abs( sum_00 );
        //corr = (abs(sum))/ (corr);
        corr = (abs(sum));
        return corr;
      }

      void print_coordinate(scitbx::af::const_ref <int> indx)
      {
        int atom_i, n_atom = indx.size();
        for(int i=0;i<n_atom;i++)
        {
          atom_i=indx[i];
          std::cout<<rtp_[atom_i][0]<<"\t"<<rtp_[atom_i][1]<<"\t"<<rtp_[atom_i][2]<<std::endl;
        }
        return;
      }

      scitbx::vec3 <FloatType> polar(scitbx::vec3<FloatType> xyz, int i)
      {
        x=xyz[0];
        y=xyz[1];
        z=xyz[2];
        r=xyz.length();
        if(r==0.0) {
        atom_cs_t_[i][0] = 1.0;
        atom_cs_t_[i][1] = 0.0;
        scitbx::vec3 <FloatType>  this_rtp(0,0,0);
        return this_rtp;
        }
        cos_t = z/r;
        t=acos(cos_t);
        atom_cs_t_[i][0] = cos_t;
        atom_cs_t_[i][1] = sin( t );
        if(x <0)
        {
         p=scitbx::constants::pi-asin(y/std::sqrt(x*x+y*y));}
        else
        {if(y<0)
                {p=scitbx::constants::pi*2.0+asin(y/std::sqrt(x*x+y*y));}
                else {p=asin(y/std::sqrt(x*x+y*y));}
        }
        scitbx::vec3 <FloatType> this_rtp(r,t,p);

        return this_rtp;
      }

    std::complex<FloatType> A_q( int q_index, FloatType theta, FloatType phi )
    {
        std::complex<FloatType> result(0,0);
        int count = 0;
        for(int l=0;l<=max_L_; l+=1)
          for(int m=-l; m<=l; m++)
          {
            result += A_array_[q_index][count]*spherical_harmonic(l,m,theta, phi);
            count++;
          }
        return result;
    }

    scitbx::af::shared<FloatType> I_q_omega_default( int q_index )
    {
        return I_q_omega( q_index, grid_.const_ref() );
    }

    scitbx::af::shared<FloatType> I_q_omega( int q_index, scitbx::af::const_ref< scitbx::vec2<FloatType> > grid)
    {
        int n = grid.size();
        scitbx::af::shared<FloatType> I_q;
        for(int i=0;i<n;i++)
        {
        std::complex<FloatType> result(0,0);
        int count = 0;
        for(int l=0;l<=max_L_; l+=1)
          for(int m=-l; m<=l; m++)
          {
            result += A_array_[q_index][count]*spherical_harmonic(l,m,grid[i][0],grid[i][1]);
            count++;
          }
        I_q.push_back( norm(result) );
        }
        return I_q;
    }

    scitbx::af::shared< scitbx::vec3<FloatType> > get_grid_xyz( FloatType q_value )
    {
        scitbx::af::shared< scitbx::vec3<FloatType> > xyz;
        for(int i=0;i<grid_.size();i++)
          xyz.push_back( cartesian( grid_[i] ) * q_value );
        return xyz;
    }

    scitbx::vec3<FloatType> cartesian(scitbx::vec2<FloatType> polar)
    {
        scitbx::vec3<FloatType> this_xyz;
        this_xyz[0] = sin(polar[0]) * cos(polar[1]);
        this_xyz[1] = sin(polar[0]) * sin(polar[1]);
        this_xyz[2] = cos(polar[0]);
        return this_xyz;
    }

    scitbx::af::shared< FloatType > simulate_correlation( int q_index, int n_phi)
    {
      int n_grid = grid_.size();
      FloatType two_pi = scitbx::constants::pi *2.0;
      FloatType phi_step = two_pi / FloatType( n_phi );
      scitbx::af::shared< FloatType > I_tp;
      scitbx::af::shared< FloatType > correlation;
      scitbx::af::shared< int > phi_array_count;
      scitbx::af::shared< FloatType > I_phi_array;
      scitbx::vec2<FloatType> cs_i, cs_j;
      FloatType p0, p1, angle, sign_of_angle;
      int indx;
      FloatType d_sqr, qq;
      qq = scat_lib_.q(q_index);  //q
      d_sqr = qq/two_pi/2.0;  //d = 0.5*q*q/kapa, but we need d^2/q^2, so one q is cancelled.
      d_sqr *= d_sqr;

      for(int i=0; i<n_phi; i++)
      {
        I_phi_array.push_back(0);
        phi_array_count.push_back(1);
      }

      for(int i=0; i < n_grid; i++)
      {
        I_tp.push_back( norm(A_q( q_index, grid_[i][0], grid_[i][1]) ) );
      }

       scitbx::af::shared< scitbx::vec2<FloatType> > grid_cs_p_;
      for(int i=0; i<n_grid;i++)
      {
        grid_cs_p_.push_back(scitbx::vec2<FloatType>( cos(grid_[i][1]), sin(grid_[i][1]) ) );
      }
      for(int i=0; i < n_grid; i++)
      {
        cs_i = grid_cs_t_[i];
        p0 = grid_[i][1];
        for(int j=0; j <n_grid; j++)
        {
          if( i == j)
            indx = 0;
          else
          {
            p1 = grid_[j][1];
            cs_j = grid_cs_t_[j];
            cos_ad = cs_i[0]*cs_j[0]+cs_i[1]*cs_j[1]*cos(p1-p0);
            sin_ad =  grid_cs_p_[i][0] * grid_cs_p_[j][0] *(cs_i[1]-cs_j[1]);
//          sin_ad =  cos(p0)* cos(p1) *(cs_i[1]-cs_j[1]);
            if( cos_ad <2*d_sqr -1)  //not possible for this q
              continue;
            cos_ad = (d_sqr-cos_ad)/(d_sqr-1.0); //testing
            angle = acos( cos_ad );
            if( sin_ad < 0)
              angle = two_pi - angle;

            //cos_ad = cos(p0)*cs_j[1]*sin(p1)-cs_i[1]*sin(p0)*cos(p1);
            //angle = std::atan2( cos_ad, sin_ad);
            //if(angle < 0)
        //      angle += two_pi;
            indx = floor( angle/phi_step+0.5);
          }
          I_phi_array[indx] +=  ( I_tp[i]*I_tp[j] );
          phi_array_count[indx] += 1;
        }
      }

      I_phi_array[0] /= FloatType( phi_array_count[0] );
      for(int i=1; i< n_phi; i++)
      {
        I_phi_array[i] /= FloatType( phi_array_count[i] );
        I_phi_array[i] /= I_phi_array[0];
        std::cout<<phi_array_count[i]<<std::endl;
      }
      I_phi_array[0] = 1;
      return I_phi_array;
    }


    private:
      model <FloatType> model_;
      scitbx::af::shared< scitbx::vec3<FloatType> > xyz_;
      scitbx::vec3<FloatType> center_;
      scitbx::af::shared< scitbx::vec3<FloatType> > rtp_;
      scitbx::af::shared< scitbx::vec2<FloatType> > grid_;
      scitbx::af::shared< scitbx::vec2<FloatType> > grid_cs_t_;
      scitbx::af::shared< scitbx::vec2<FloatType> > atom_cs_t_;
      scitbx::af::shared <FloatType> F_array_;
      scitbx::af::shared< FloatType > radius_;
      scattering_library< FloatType > scat_lib_;
      scitbx::af::shared< integrated <FloatType> > Integral_;
      scitbx::af::shared< FloatType > I_all_;
      scitbx::af::shared< FloatType > I_A_;
      scitbx::af::shared< FloatType > I_B_;
      scitbx::af::shared< FloatType > I_C_;
      int n_, max_i_, max_L_;
      FloatType max_z_;
      FloatType f_w_max_, f_w_min_,q_step_,f_resl_;
      FloatType upper_limit_, lower_limit_;
      FloatType I_, r_w, delta_, rho_, d_rho_;
      //besselr<FloatType> bessel_array_;
      scitbx::af::shared< scitbx::vec3<int > > neighbors_;
      scitbx::af::shared< scitbx::af::shared <std::complex<FloatType> > > A_array_;
      scitbx::af::shared< scitbx::af::shared <std::complex<FloatType> > > C_array_;
      scitbx::af::shared< scitbx::af::shared <std::complex<FloatType> > > B_array_;
      scitbx::af::shared< scitbx::af::shared <std::complex<FloatType> > > S_array_;
      scitbx::af::shared< scitbx::af::shared <std::complex<FloatType> > > I_lm_array_;
      scitbx::af::shared< FloatType > wigner_array_;
      scitbx::af::shared< FloatType > coef_array_;
      bool add_wigner;
      int accum_, coef_accum_;
      //scitbx::math::acos_table<FloatType> acos_table_;
      FloatType distance, cos_ad, sin_ad, r1,t1,p1,t0,p0, ct0, st0;
      FloatType r,t,p, x,y,z, cos_t;
      scitbx::af::shared< scitbx::af::shared< scitbx::af::shared< FloatType > > > spatial_correlation_;
      scitbx::af::shared< scitbx::af::shared< FloatType > > expansion_coef_array_;
      wigner3j_fast<FloatType> wigner_fast_;
  };




}}  // namespace sastbx::intensity
#endif // SASTBX_SHE_H
