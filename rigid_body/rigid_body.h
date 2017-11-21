#ifndef SASTBX_RB_H
#define SASTBX_RB_H
#endif

#include <scitbx/constants.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>

#include <scitbx/array_family/shared.h>
#include <map>

#include <string>
#include <iomanip>
#include <complex>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <scitbx/mat3.h>
#include <cctbx/sgtbx/rt_mx.h>
#include <scitbx/histogram.h>
#include <scitbx/math/r3_rotation.h>

using namespace std;
namespace sastbx { namespace rigid_body {

   template<typename FloatType>
   class rigidbody 
   {
        public:
	  rigidbody( scitbx::af::const_ref< scitbx::vec3< FloatType> > const& coords,
		     scitbx::af::const_ref< int > const& sample_indices,
		     FloatType dmax,
		     int max_i,
		     bool griding
		   ): dmax_(dmax), angle_dev_(3.0), center_(0,0,0), max_i_(max_i), griding_(griding), buffer_(0.0), layer_thick_(1.8), w1_(1.0), w2_(1.0), w3_(1.0), dist_cutoff_(10.0), dist_cutoff_sq_(dist_cutoff_*dist_cutoff_)
          {
		natom_ = coords.size();
		nsample_ = sample_indices.size();
		for(int i=0;i<natom_;i++)
		  coords_.push_back(coords[i]); 
 
		for(int i=0;i<nsample_;i++)
		{
		   sample_coords_.push_back( coords[ sample_indices[i] ] );
		   center_ += sample_coords_[ i ];
		}
		center_ = center_ /nsample_;
		//shift to the center_
		FloatType factor = sqrt( 2.0 * 3.1415926 );
		angle_dev_ = angle_dev_ * 3.1415926 / 180.0;
		for(int i=0;i<nsample_;i++)
		{
		  sample_coords_[i] -= center_;
		}
		for(int i=0; i< natom_; i++)
		  coords_[i] -= center_;

		new_coords_=sample_coords_.deep_copy();
		backup_coords_ = sample_coords_.deep_copy();
        //        if(griding_)		
                  build_grid();
		pick_surface_atoms(buffer_);
	  }

	  void reset()
	  {
		sample_coords_.clear();
		sample_coords_ = backup_coords_.deep_copy();
	  }

	  scitbx::af::shared< scitbx::vec3< FloatType > > rotate_translate( scitbx::vec3< FloatType > vector, scitbx::vec3< FloatType > angle)
	  {
		scitbx::mat3< FloatType > rotation_matrix = euler_xyz_matrix( angle );
		//shift_ = center_ + vector;
		shift_ = vector;
                for(int i=0;i<nsample_;i++)
                  new_coords_[i] = rotation_matrix * sample_coords_[i] + shift_;
		return new_coords_;
	  } 


          scitbx::af::shared< scitbx::vec3< FloatType > > rotate_translate_q( scitbx::vec3< FloatType > vector, scitbx::af::const_ref<FloatType> q)
          {
                scitbx::mat3< FloatType > rotation_matrix = quaternion_matrix(q);
		shift_ = vector;
                //shift_ = center_ + vector;
                for(int i=0;i<nsample_;i++)
                  new_coords_[i] = rotation_matrix * sample_coords_[i] + shift_;
                return new_coords_;
          }

          scitbx::af::shared< scitbx::vec3< FloatType > > rotate_translate_v( scitbx::vec3< FloatType > vector,scitbx::vec3< FloatType > vector_r, FloatType angle)
          {
                scitbx::mat3< FloatType > rotation_matrix = scitbx::math::r3_rotation::axis_and_angle_as_matrix(vector_r, angle);
                //shift_ = center_ + vector;
		shift_ = vector;
                for(int i=0;i<nsample_;i++)
                  new_coords_[i] = rotation_matrix * sample_coords_[i] + shift_;
                return new_coords_;
          }

	  void translate_only( scitbx::vec3< FloatType > vector )
	  {
	//	shift_ = center_ + vector;
		shift_ = vector;
		for(int i=0; i< nsample_; i++)
		  new_coords_[i] = sample_coords_[i] + shift_;
	  }

	  void rotate_only( scitbx::vec3< FloatType > vector, FloatType angle )
	  {
		scitbx::mat3< FloatType > rotation_matrix = scitbx::math::r3_rotation::axis_and_angle_as_matrix(vector, angle);
	        for( int i = 0; i < nsample_; i++)
		  sample_coords_[i] = rotation_matrix * sample_coords_[i];
	  }

          /* rotation around given anchor point
           * The anchor point coordinate should be relative to the center of the corresponding rigid-body
           * in other words, point = original_point-center_
           */
          scitbx::af::shared< scitbx::vec3< FloatType > > rotate_around_one_point( scitbx::vec3< FloatType > point, scitbx::vec3< FloatType > angles )
          {
                scitbx::mat3< FloatType > rotation_matrix = euler_xyz_matrix( angles );
                shift_ = point;
                for(int i=0;i<nsample_;i++)
                  new_coords_[i] = rotation_matrix * (sample_coords_[i] + shift_);
                return new_coords_;
          }


          /* rotation around given two anchor points
           * The reference of the anchor points are the center_
           * the relative location of the two points set a vector
           * the rotation will be around that vector
           */
          scitbx::af::shared< scitbx::vec3< FloatType > > rotate_around_two_point( scitbx::vec3< FloatType > point1, scitbx::vec3< FloatType > point2,FloatType angle )
          {
	        scitbx::vec3< FloatType > vector;
                vector = point2-point1;
                scitbx::mat3< FloatType > rotation_matrix = scitbx::math::r3_rotation::axis_and_angle_as_matrix(vector, angle);
                for(int i=0;i<nsample_;i++)
                  new_coords_[i] = rotation_matrix * (sample_coords_[i] + point1);
                return new_coords_;
          }

	  scitbx::mat3<FloatType> quaternion_matrix( scitbx::af::const_ref< FloatType > q) 
	  {
	    return scitbx::math::r3_rotation::unit_quaternion_as_matrix(q[0],q[1],q[2],q[3]);
	  }

	  scitbx::mat3<FloatType> euler_xyz_matrix( scitbx::vec3< FloatType > ea)
	  {
            // The euler angles are in unit of radiant
		FloatType cx,sx,cy,sy,cz,sz;

		cx = cos(ea[0]); 
		sx = sin(ea[0]);
  		cy = cos(ea[1]);
  		sy = sin(ea[1]);
  		cz = cos(ea[2]);
  		sz = sin(ea[2]);

	  return scitbx::mat3< FloatType> (
		  cy*cz,         -cy*sz,     sy,
		  cz*sx*sy+cx*sz, cx*cz-sx*sy*sz, -cy*sx,
    		  -cx*cz*sy+sx*sz, cz*sx+cx*sy*sz,  cx*cy);

	  }	  

  
          scitbx::af::shared < scitbx::vec3<FloatType> > get_crd()
  	  {
		   return new_coords_;
	  }
          void calc_histogram()
          {
                int count=0, surface_i, surface_j, indx_i, indx_j;
		FloatType double_layer;
		double_layer = layer_thick_*2.0;
                scitbx::vec3 <FloatType> Vdiff;
                FloatType dist;
                indx_i = 0;
                for(int i=0;i<nsample_-1;i++)
                {
                  for(int j=i+1;j<nsample_;j++)
                   {
                        Vdiff = sample_coords_[i]-sample_coords_[j];
                        dist = sqrt( Vdiff*Vdiff );
                        dist_array_p2p_.push_back( dist );
                        if(surface_[i]) {
                          dist_array_p2s_.push_back( dist+layer_thick_ );
                          if(surface_[j])
                          { dist_array_s2s_.push_back( dist+double_layer );
                            //dist_array_.push_back( dist ); //this does not improve the fitting
                            dist_array_p2s_.push_back( dist+layer_thick_ );
                          }
                        }
                        else if(surface_[j]) { dist_array_p2s_.push_back( dist+layer_thick_ );}
                   }
                   if(surface_[i]) dist_array_p2s_.push_back( layer_thick_ );
                }

                int n_slots=static_cast<int> (dmax_ + 0.5);
                scitbx::histogram<FloatType, int> h1(dist_array_p2p_.const_ref(),0.0, dmax_, n_slots);
                raw_h1_ = h1.slots();
                scitbx::histogram<FloatType, int> h2(dist_array_p2s_.const_ref(),0.0, dmax_, n_slots);
                raw_h2_ = h2.slots();
                scitbx::histogram<FloatType, int> h3(dist_array_s2s_.const_ref(),0.0, dmax_, n_slots);
                raw_h3_ = h3.slots();

                histogram_p2p_.push_back(raw_h1_[0]*0.75 + raw_h1_[1]/4.0 );
                histogram_p2s_.push_back(raw_h2_[0]*0.75 + raw_h2_[1]/4.0 );
                histogram_s2s_.push_back(raw_h3_[0]*0.75 + raw_h3_[1]/4.0 );
                for(int i=1;i<n_slots-1;i++){
                  histogram_p2p_.push_back((raw_h1_[i-1]+raw_h1_[i+1])/4.0 + raw_h1_[i]/2.0);
                  histogram_p2s_.push_back((raw_h2_[i-1]+raw_h2_[i+1])/4.0 + raw_h2_[i]/2.0);
                  histogram_s2s_.push_back((raw_h3_[i-1]+raw_h3_[i+1])/4.0 + raw_h3_[i]/2.0);
		}
                histogram_p2p_.push_back(raw_h1_[n_slots-1]*0.75+raw_h1_[n_slots-2]/4.0);
                histogram_p2s_.push_back(raw_h2_[n_slots-1]*0.75+raw_h2_[n_slots-2]/4.0);
                histogram_s2s_.push_back(raw_h3_[n_slots-1]*0.75+raw_h3_[n_slots-2]/4.0);

          }
/* 
	  void calc_histogram()
	  {
          	int count=0, surface_i, surface_j, indx_i, indx_j;
          	scitbx::vec3 <FloatType> Vdiff;
		FloatType dist;
		indx_i = 0;
          	for(int i=0;i<nsample_-1;i++)
          	{ 
		  for(int j=i+1;j<nsample_;j++)
            	   {
                	Vdiff = sample_coords_[i]-sample_coords_[j];
			dist = sqrt( Vdiff*Vdiff );
                	dist_array_.push_back( dist );
            	   }
          	}
	
		int n_water = dummy_wat_.size();
		std::cout<< n_water << std::endl;

                for(int i=0;i<n_water-1;i++)
                {
                  for(int j=i+1;j<n_water;j++)
                   {
                        Vdiff = dummy_wat_[i]-dummy_wat_[j];
                        dist = sqrt( Vdiff*Vdiff );
                        dist_array_.push_back( dist );
                   }
                }

                for(int i=0;i<n_water;i++)
                {
                  for(int j=0;j<nsample_;j++)
                   {
                        Vdiff = new_coords_[i]-dummy_wat_[j];
                        dist = sqrt( Vdiff*Vdiff );
                        dist_array_.push_back( dist );
                   }
                }


		int n_slots=static_cast<int> (dmax_ + 0.5);
		scitbx::histogram<FloatType, int> h(dist_array_.const_ref(),0.0, dmax_, n_slots);		
		raw_h_ = h.slots();
		histogram_.push_back(raw_h_[0]*0.75 + raw_h_[1]/4.0 );
		for(int i=1;i<n_slots-1;i++)
	          histogram_.push_back((raw_h_[i-1]+raw_h_[i+1])/4.0 + raw_h_[i]/2.0);
        	histogram_.push_back(raw_h_[n_slots-1]*0.75+raw_h_[n_slots-2]/4.0);

	  }
*/
	  scitbx::af::shared< FloatType >  that_histogram( scitbx::af::shared< scitbx::vec3< FloatType > > that_crd)
	  {
          	int count=0;
		int that_size = that_crd.size();
          	scitbx::vec3 <FloatType> Vdiff;
		int n_slots=static_cast<int> (dmax_ + 0.5);

		scitbx::af::shared< FloatType > tmp_array(nsample_*that_size, scitbx::af::init_functor_null<FloatType>() );

          	for(int i=0;i<nsample_;i++)
          	{  for(int j=0;j<that_size;j++)
            	   {
                	Vdiff = new_coords_[i]-that_crd[j];
                	tmp_array[count] = sqrt(Vdiff*Vdiff);
                	count++;
            	   }
          	}
		

		scitbx::af::shared< FloatType > that_hist; //n_slots, scitbx::af::init_functor_null<FloatType>() );

		scitbx::histogram<FloatType, int> h(tmp_array.const_ref(),0.0, dmax_, n_slots);	
		raw_h_ = h.slots();
		that_hist.push_back(raw_h_[0]*0.75 + raw_h_[1]/4.0 );
		for(int i=1;i<n_slots-1;i++)
	          that_hist.push_back((raw_h_[i-1]+raw_h_[i+1])/4.0 + raw_h_[i]/2.0);
        	that_hist.push_back(raw_h_[n_slots-1]*0.75+raw_h_[n_slots-2]/4.0);

            return that_hist;

	  }
	  
	  scitbx::af::shared< FloatType >  get_hist_var()
	  {

	   return hist_var_;

	  }

          int get_clash()
	  { return raw_h1_[0] + raw_h1_[1];
	  }

	  void set_weights(FloatType w2, FloatType w3)
	  {
		w2_=w2;
		w3_=w3;
		return;
	  }

	  scitbx::af::shared<FloatType> get_histogram()
	  {
		if(histogram_.size() == 0)
		  calc_histogram();
		histogram_.clear();
		for(int i=0; i< histogram_p2p_.size(); i++)
		{
		 histogram_.push_back(histogram_p2p_[i] + histogram_p2s_[i]*w2_ + histogram_s2s_[i]*w3_);
	        }
		return histogram_;
	  }

      void build_grid()
      {
        int i,max_grid,second_max;
        FloatType theta, phi, f_max, two_pi;

        max_grid=fibonacci(max_i_);
        second_max=fibonacci(max_i_-1);
        f_max = static_cast<FloatType> (max_grid);
        two_pi = scitbx::constants::pi * 2.0;

        for(i=0;i<natom_;i++)
	  rtp_.push_back( polar( coords_[i] ) );

        for(i=0;i<=max_grid;i++)
         {
                scitbx::vec3<FloatType> omega;
                theta =acos(1.0-2.0*i/f_max); //[0,pi]
                phi = two_pi * ((i*max_grid) % second_max)/static_cast<FloatType>(second_max);
                omega[0]=theta;
                omega[1]=phi;
		omega[2]=0;
                grid_.push_back(omega);
         }
	F_omega_array();
        return ;
      }//end of Griding

      
      void F_omega_array() // find the surface vector, as described in CRYSOL paper
      {
        int i, total, max_grid, indx, second_max,span;

        total=grid_.size();

        max_grid=total-1;
        second_max=fibonacci(max_i_-1);
        span=max_i_*4.0;

        for(i=0;i<natom_;i++)
        {
          for(int j=0;j<total;j++)
          {
                update_F_w(j,i);
          }
        }
      // Fill some gaps, in case no atom is found in that direction
        for(int j=1;j<total-1;j++)
        { if(grid_[j][2]==0.0) grid_[j][2]=(grid_[j-1][2]+grid_[j+1][2])*0.5;}
        if(grid_[0][2]==0.0) grid_[0][2]=grid_[1][2];
        if(grid_[total-1][2]==0.0) grid_[total-1][2]=grid_[total-2][2];

	f_w_min_ = 100000;
	f_w_max_ = -1.0;
        for(int j=0;j<total;j++) {
          if(f_w_min_> grid_[j][2]) f_w_min_=grid_[j][2];
          if(f_w_max_< grid_[j][2]) f_w_max_=grid_[j][2];
        }
        return ;
      }

      void update_F_w(int direction_j, int atom_i)
      {
        FloatType distance, cos_ad, sin_ad;
        FloatType r1,t1,p1,t0,p0;
        t0=grid_[direction_j][0];
        p0=grid_[direction_j][1];

        r1=rtp_[atom_i][0];
        t1=rtp_[atom_i][1];
        p1=rtp_[atom_i][2];
        cos_ad = cos(t0)*cos(t1)+sin(t0)*sin(t1)*cos(p1-p0);
        if(cos_ad >=0)
        {
          sin_ad = sqrt(1.0 - cos_ad*cos_ad);
          distance = r1*sin_ad;
          if(distance < 3.5 )  //+radius_[atom_i])) //need to be corrected to r_w + r_[i]
          {
             distance = r1*cos_ad;
             if(grid_[direction_j][2] < distance+ 1.0)
             { grid_[direction_j][2] = distance + 1.0; }
          }
        }
        return;
      }

      scitbx::af::shared< scitbx::vec3< FloatType> > get_grid()
      {
	return grid_;
      }

      FloatType range_f_w()
      {
	return f_w_max_; // - f_w_min_;
      }
      
      scitbx::vec3 <FloatType> polar(scitbx::vec3<FloatType> xyz)
      {
        FloatType r,t,p;
        FloatType x,y,z;
        x=xyz[0];
        y=xyz[1];
        z=xyz[2];
        r=xyz.length();
        if(r==0.0) {
        scitbx::vec3 <FloatType>  this_rtp(0,0,0);
        return this_rtp;
        }

        t=acos(z/r);
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

      scitbx::vec3< FloatType > get_center()
      { return center_; }

      int fibonacci(int n)
      {
        if(n==0) {return 0;}
        else if(n == 1) { return 1;}
        else { return fibonacci(n-1) + fibonacci(n-2);}
      }

     scitbx::af::shared< int > get_surface_atoms()
     {
	return surface_atoms_;
     }

     void pick_surface_atoms(FloatType cutoff)
     {
	xmin_ = xmax_ = coords_[0][0];
	ymin_ = ymax_ = coords_[0][1];
	zmin_ = zmax_ = coords_[0][2];
	for( int i=1; i< natom_; i++)
	{
	  if( coords_[i][0] > xmax_ ) xmax_=coords_[i][0];
	  if( coords_[i][1] > ymax_ ) ymax_=coords_[i][1];
	  if( coords_[i][2] > zmax_ ) zmax_=coords_[i][2];
	  if( coords_[i][0] < xmin_ ) xmin_=coords_[i][0];
	  if( coords_[i][1] < ymin_ ) ymin_=coords_[i][1];
	  if( coords_[i][2] < zmin_ ) zmin_=coords_[i][2];
	} 
	int nx,ny,nz, nyz, nxyz;
	int buffer;
	int this_x, this_y, this_z, indx;
 	int new_indx, i,j,k, ii, jj,kk;
        FloatType size;
        size = 3.0;

	buffer = 10;
	nx = int((xmax_-xmin_+0.5)/size) + buffer;
	ny = int((ymax_-ymin_+0.5)/size) + buffer;
	nz = int((zmax_-zmin_+0.5)/size) + buffer;
	buffer = buffer / 2;
	scitbx::af::shared < int  > flag( nx*ny*nz, 0 );
	scitbx::af::shared< scitbx::vec3< FloatType > > atom_indices;
	nyz = ny*nz;
	nxyz = nx*nyz;
	scitbx::af::shared<FloatType> shifts ;
	shifts.push_back( nyz );
	shifts.push_back( -nyz );
	shifts.push_back( nz );
	shifts.push_back( -nz );
	shifts.push_back( 1 );
	shifts.push_back( -1 );

	for( i =0; i< natom_; i++)
	{
	  this_x = int((coords_[i][0]-xmin_)/size) + buffer;
	  this_y = int((coords_[i][1]-ymin_)/size) + buffer;
	  this_z = int((coords_[i][2]-zmin_)/size) + buffer;
	  atom_indices.push_back(scitbx::vec3<FloatType> (this_x, this_y, this_z) );

          indx = this_x *nyz + this_y * nz + this_z;
	  for(ii = 0; ii<6; ii++)
	  { new_indx = indx + shifts[ii];
          flag[ new_indx ] = 1;
	  }
	  
/*	  for( ii =-1; ii<2; ii++)
	  {  for( jj = -1; jj < 2; jj++)
           {   for(kk = -1; kk < 2; kk++)
	      {
		indx = (this_x + ii) *nyz + (this_y + jj) * nz + (this_z + kk);
	        flag[ indx ] = 1;
		//std::cout<<this_x<<" "<<this_y<<" "<<this_z<<std::endl;
	      }
	   }
	   }
*/	}

        bool marked;
	scitbx::vec3<FloatType> min_xyz(xmin_, ymin_, zmin_);

	for( i = 2 ; i < nx-2; i++)
	for( j = 2 ; j < ny-2; j++)
	for( k = 2 ; k < nz-2; k++)
	{
	  indx = i * nyz + j * nz + k;
	  marked = false;
	  if(flag[indx]==1)
	  {
	   for( ii = 0; ii < 6; ii++)
	   {
		new_indx = indx + shifts[ii];
                if(flag[ new_indx ] == 0)
                {  flag[indx]=2;
                  scitbx::vec3<FloatType> tmp_xyz( (i-buffer)*size,(j-buffer)*size, (k-buffer)*size);
                  tmp_xyz = tmp_xyz + min_xyz;
                  dummy_wat_.push_back( tmp_xyz );
	 	}
	   }
/*	  for( ii =-1; ii<2; ii++)
            {for( jj = -1; jj < 2; jj++)
              {for(kk = -1; kk < 2; kk++)
              {
                new_indx = indx + ii *nyz + jj * nz + kk;
                if(flag[ new_indx ] == 0)
		{  flag[indx]=2;
		  scitbx::vec3<FloatType> tmp_xyz( (i-buffer)*size,(j-buffer)*size, (k-buffer)*size);
		  tmp_xyz = tmp_xyz + min_xyz;
		  dummy_wat_.push_back( tmp_xyz );
		  marked = true;
		}
	      if(marked) break;
              }
	      if(marked) break;
              }
	     if(marked) break;
	     }//end of ii
*/
/*
		std::cout<<"ATOM   1400  N   ILE A 178    ";
		std::cout<<setw(8)<<setprecision(3)<<(i-buffer)*size+xmin_+center_[0];
		std::cout<<setw(8)<<setprecision(3)<<(j-buffer)*size+ymin_+center_[1];
		std::cout<<setw(8)<<setprecision(3)<<(k-buffer)*size+zmin_+center_[2];scitbx::histogram<FloatType, int> h(tmp_array.const_ref(),0.0, dmax_, n_slots);
		std::cout<<setw(5)<<setprecision(2)<<1.0;
		std::cout<<setw(5)<<setprecision(2)<<flag[indx];
		std::cout<<std::endl;
*/	   } //endif
	}
 
	bool surface;
	for( int i=0; i < natom_; i++)
	{
	  this_x = atom_indices[i][0];
	  this_y = atom_indices[i][1];
	  this_z = atom_indices[i][2];
	  indx = this_x*nyz + this_y*nz + this_z;
	  surface = false;
	  for( ii = 0; ii < 6; ii++)
	  {
		new_indx = indx + shifts[ii];
		if( flag[new_indx] == 2)
		{
		  surface = true;
		  surface_atoms_.push_back(i);
		  break;
		}
	  }
/*          for( ii =-1; ii<2; ii+=2)
	  {
            for( jj = -1; jj < 2; jj+=2)
            {  for(kk = -1; kk < 2; kk+=2)
              {
                //indx = (this_x + ii*2) *nyz + (this_y + jj*2) * nz + (this_z + kk*2);
                new_indx = indx + ii*nyz + jj*nz + kk;
                if( flag[ new_indx ] == 2 )
		{
		  surface = true;
                  surface_atoms_.push_back(i);
		  break;
		}
              }
	      if(surface) break;
	      }
	    if(surface) break;
	   }
*/
	  surface_.push_back( surface );

	}
	return ;
     }

     scitbx::af::shared< scitbx::vec2< FloatType > > update_contact_list( scitbx::af::const_ref< scitbx::vec3<FloatType> > that_crd )
     {
	int that_size, i,j;
	scitbx::vec3<FloatType> Vdiff;
	FloatType tmp_distsq;
	contact_list_.clear();
	that_size = that_crd.size();	
	for(i=0; i < natom_; i ++)
	  for(j = 0; j<that_size; j++)
	  {
	     Vdiff = new_coords_[i]-that_crd[j];
	     tmp_distsq = Vdiff*Vdiff;
	     if( tmp_distsq < dist_cutoff_sq_)
	     {
		contact_list_.push_back( scitbx::vec2<FloatType> (i,j) );
		//dist_array_contact_.push_back( sqrt(tmp_distsq) );
	     }
	  }
	return contact_list_;
     }

     scitbx::af::shared<int>  calc_contact_hist(scitbx::af::const_ref< scitbx::vec3<FloatType> > that_crd, scitbx::af::const_ref< scitbx::vec2<FloatType > > contact_list)
     {
	int list_size,i,j,k;
	scitbx::vec3<FloatType> Vdiff;
	scitbx::vec2<FloatType> pair;
	list_size = contact_list.size();
	dist_array_contact_.clear();
	
	if(list_size == 0) update_contact_list( that_crd );
	
	for(k =0;k<list_size; k++)
	{
	   pair = contact_list[k];
	   Vdiff = new_coords_[ pair[0] ] - that_crd[ pair[1] ];
	   dist_array_contact_.push_back( sqrt( Vdiff*Vdiff ) );
	}
	scitbx::histogram<FloatType, int> h(dist_array_contact_.const_ref(),0.0, 5.0, 5);
	
	return h.slots();
	
     }
     
     


          
        private:
	  FloatType xmin_, ymin_, zmin_, xmax_, ymax_, zmax_;
	  scitbx::af::shared< scitbx::vec3< FloatType > > coords_;
	  scitbx::af::shared< scitbx::vec3< FloatType > > sample_coords_;
	  scitbx::af::shared< scitbx::vec3< FloatType > > backup_coords_;
          scitbx::af::shared< scitbx::vec3< FloatType > > new_coords_; 
	  scitbx::af::shared< FloatType > r_;
	  scitbx::af::shared< FloatType > sigma_;
	  scitbx::af::shared< FloatType > p_;
	  scitbx::af::shared< FloatType > hist_var_;
	  scitbx::af::shared< FloatType > dist_array_p2p_;
	  scitbx::af::shared< FloatType > dist_array_p2s_;
	  scitbx::af::shared< FloatType > dist_array_s2s_;
	  scitbx::af::shared< FloatType > histogram_;
	  scitbx::af::shared< FloatType > histogram_p2p_;
	  scitbx::af::shared< FloatType > histogram_p2s_;
	  scitbx::af::shared< FloatType > histogram_s2s_;
	  scitbx::af::shared< scitbx::vec3< FloatType > > rtp_;
          scitbx::af::shared< scitbx::vec3< FloatType > > grid_; 
          scitbx::af::shared< scitbx::vec3< FloatType > > dummy_wat_; 
          scitbx::af::shared< scitbx::vec2< FloatType > > contact_list_;
	  scitbx::af::shared< FloatType > dist_array_contact_;
	  scitbx::af::shared< int > raw_h1_, raw_h2_, raw_h3_, raw_h_;
	  scitbx::af::shared< int > surface_atoms_;
	  scitbx::af::shared< bool > surface_;
	  scitbx::vec3< FloatType > center_, shift_;
	  FloatType w1_, w2_, w3_;
	  int natom_, nsample_, max_i_;
	  FloatType dmax_, angle_dev_, f_w_max_, f_w_min_, buffer_, layer_thick_, dist_cutoff_, dist_cutoff_sq_;
	  bool griding_;

   }; //end of rigidbody class



} //namepsace rigid_body
} //namespace sastbx

