#include <sastbx/fXS/image_simulator_tools.h>
#if __APPLE__
void sincos(double t, double *s, double*c) {*s=sin(t); *c=cos(t);}
#endif

namespace sastbx {
namespace fXS {

  /* ==========================================================================
     image_composer Class

     This class pre-calculates the Miller indices required to construct a
     scattering pattern with the given parameters.  This class is used in
     conjunction with the image_simulator class in Python.

     Public Methods:
       image_composer - 2 constructors
       set_detector_size - self-explanatory (scitbx::vec2<int>)
       set_beam_center - self-explanatory (scitbx::vec2<int>)
       set_pixel_size - self-explanatory
       set_resolution - self-explanatory
       get_resolution - self-explanatory
       set_ewald_sphere - self-explanatory (sastbx::fXS::ewald_sphere)
       cache_hkl - caches the Miller index of each pixel and returns an array
                   of integral Miller indices required for interpolation
       build_image - given an array of intensities, returns an array of pixel
                     intensities scaled to a maximum value
       get_q - returns the q value for each pixel

     --------------------------------------------------------------------------
     Constructor

     Arguments:
       None

       ds    - detector size (scitbx::vec2<int>)
       bc    - beam center (scitbx::vec2<int>)
       ps    - pixel size (double)
       es    - ewald sphere (sastbx::fXS::ewald_sphere)

     Returns:
       None

     Notes:
       Every method in the non-trivial constructor needs to be called before
         building an image.
     --------------------------------------------------------------------------
  */
  sastbx::fXS::image_composer::image_composer() {}

  sastbx::fXS::image_composer::image_composer
  (const scitbx::vec2<int>& ds,const scitbx::vec2<int>& bc,
   const double& ps,const sastbx::fXS::ewald_sphere& es) {
    set_detector_size(ds);
    set_beam_center(bc);
    set_pixel_size(ps);
    set_ewald_sphere(es);
  }

  void sastbx::fXS::image_composer::set_detector_size
  (const scitbx::vec2<int>& ds) {
    detector_size = ds;
  }

  scitbx::vec2<int> sastbx::fXS::image_composer::get_detector_size() {
    return detector_size;
  }

  void sastbx::fXS::image_composer::set_beam_center
  (const scitbx::vec2<int>& bc) {
    beam_center = bc;
  }

  scitbx::vec2<int> sastbx::fXS::image_composer::get_beam_center() {
    return beam_center;
  }

  void sastbx::fXS::image_composer::set_pixel_size(const double& ps) {
    pixel_size = ps;
  }

  double sastbx::fXS::image_composer::get_pixel_size() {
    return pixel_size;
  }

  void sastbx::fXS::image_composer::set_ewald_sphere
  (const sastbx::fXS::ewald_sphere& es) {
    ewald_sphere = es;
  }

  sastbx::fXS::ewald_sphere sastbx::fXS::image_composer::get_ewald_sphere() {
    return ewald_sphere;
  }

  double sastbx::fXS::image_composer::get_distance() {
    return ewald_sphere.get_distance();
  }

  /* --------------------------------------------------------------------------
     Caches the Miller indices required for generating scattering pattern

     Arguments:
       None

     Returns:
       array of integer Miller indices (scitbx::af::shared
                                        <scitbx::vec3<double> >)

     Notes:
       This method sets the private data member, h_cache.
       h_cache stores the matching reciprocal space vector for every pixel
     --------------------------------------------------------------------------
  */
  scitbx::af::shared<scitbx::vec3<double> >
  sastbx::fXS::image_composer::cache_h() {
    int n_x = detector_size[0] + 1;
    int n_y = detector_size[1] + 1;
    int n_corners = n_x * n_y;
    h_cache.clear();
    h_cache.reserve(n_corners);
    scitbx::vec3<double> h;
    scitbx::vec2<double> xy;
    for (int y=0; y<n_y; y++) {
      for (int x=0; x<n_x; x++) {
        // find corners of pixels
        xy[0] = (x - beam_center[0]) * pixel_size;
        xy[1] = (y - beam_center[1]) * pixel_size;
        h = ewald_sphere.get_h(xy);
        h_cache.push_back(h);
      }
    }
    return h_cache;
  }

  /* --------------------------------------------------------------------------
     Converts structure factor intensities into pixel intensities

     Arguments:
       intensities - array of structure factor intensities with integer Miller
                     indices (scitbx::af::const_ref)

     Returns:
       array of pixel intensities (scitbx::af::shared<double>)

     Notes:
       None

     --------------------------------------------------------------------------
  */
  scitbx::af::shared<double>
  sastbx::fXS::image_composer::build_image
  (const scitbx::af::const_ref<double>& intensities) {

    // integrate intensity over area of pixel
    assert(intensities.size() == h_cache.size());
    int dsx = detector_size[0] + 1;
    scitbx::af::shared<double> pixels(detector_size[0]*detector_size[1]);
    double I00,I01,I10,I11, pixel_area;
    pixel_area = pixel_size*pixel_size;
    for (int y=0; y<detector_size[1]; y++) {
      for (int x=0; x<detector_size[0]; x++) {
        I00 = intensities[y*dsx + x];
        I01 = intensities[y*dsx + (x+1)];
        I10 = intensities[(y+1)*dsx + x];
        I11 = intensities[(y+1)*dsx + (x+1)];
        pixels[y*detector_size[0]+x] = 0.25*(I00+I01+I10+I11)*pixel_area;
      }
    }

    return pixels;
  }

  scitbx::af::shared<double> sastbx::fXS::image_composer::get_q() {
    scitbx::af::shared<double> q;
    q.reserve(detector_size[0] * detector_size[1]);
    scitbx::vec3<double> h;
    scitbx::vec2<double> xy;
    double shift = 0.5*pixel_size;
    for (int y=0; y<detector_size[1]; y++) {
      for (int x=0; x<detector_size[0]; x++) {
        // find center of pixel for q
        xy[0] = (x - beam_center[0]) * pixel_size + shift;
        xy[1] = (y - beam_center[1]) * pixel_size + shift;
        h = ewald_sphere.get_h(xy);
        q.push_back(ewald_sphere.h_to_q(h));
      }
    }
    return q;
  }

  scitbx::af::shared<double> sastbx::fXS::image_composer::get_phi() {
    scitbx::af::shared<double> phi;
    phi.reserve(detector_size[0] * detector_size[1]);
    double shift = 0.5;     // q is for center of pixel
    double xx,yy,tmp_phi;
    for (int y=0; y<detector_size[1]; y++) {
      for (int x=0; x<detector_size[0]; x++) {
        xx = x - beam_center[0] + shift;
        yy = y - beam_center[1] + shift;
        tmp_phi = std::atan2(yy,xx);     // xx == 0 should never happen
        if (tmp_phi < 0) {
          tmp_phi += scitbx::constants::two_pi;
        }
        phi.push_back(tmp_phi);
      }
    }
    return phi;
  }

  /* ==========================================================================
     end of image_composer class
  */

  /* ==========================================================================
     mean_and_variance_by_q Function

     Calculates the mean and variance for each q

     Ideally, statistical methods from scitbx::math should be used, but the
     binned_data object is weird when constructed using
     scitbx::af::shared<scitbx::af::shared<double> >.  Every bin is the same.
     --------------------------------------------------------------------------
  */
  scitbx::af::shared<scitbx::vec2<double> > mean_and_variance_by_q
  (const double& q_min,const double& q_max,const double& dq,
   const int& n_bins,const scitbx::af::const_ref<double>& q,
   const scitbx::af::const_ref<double>& data) {

    // bin data
    std::vector<std::vector<double> > binned_data(n_bins);
    int q_i;
    for (int i=0; i<q.size(); i++) {
      if ((q[i] >= q_min) && (q[i] < q_max)) {
        q_i = int((q[i] - q_min)/dq);
        binned_data[q_i].push_back(data[i]);
      }
    }

    // calculate mean and variance of each bin
    double d,sum;
    scitbx::af::shared<scitbx::vec2<double> > mv(n_bins);
    for (int i=0; i<n_bins; i++) {
      if (binned_data[i].size() > 0) {
        sum = 0.0;
        for (int j=0; j<binned_data[i].size(); j++) {
          sum += binned_data[i][j];
        }
        mv[i][0] = sum/binned_data[i].size();
        sum = 0.0;
        for (int j=0; j<binned_data[i].size(); j++) {
          d = binned_data[i][j] - mv[i][0];
            sum += d*d;
        }
        mv[i][1] = sum/binned_data[i].size();
      }
    }

    return mv;
  }

  /* ==========================================================================
     I(q)
     --------------------------------------------------------------------------
  */
  double I_q(const scitbx::af::const_ref<double>& q,
             const scitbx::af::const_ref<double>& data,
             const double& q_i, const double& dq) {
    double mean = 0.0;
    double count = 0.0;
    double q_min = q_i - 0.5*dq;
    double q_max = q_i + 0.5*dq;
    for (int i=0; i<q.size(); i++) {
      if ((q[i] >= q_min) && (q[i] < q_max)) {
        mean += data[i];
        count += 1.0;
      }
    }
    mean = mean/count;
    return mean;
  }

  /* ==========================================================================
     A(q) = sum( f(q) exp(i q * r) )
     |h| = 2 * stol, stol = sin(theta)/lambda
     |q| = 4 pi stol

     --------------------------------------------------------------------------
  */
  scitbx::af::shared<std::complex<double> > direct_sum_structure_factors
  (const scitbx::af::const_ref<std::string>& scatterers,
   const scitbx::af::const_ref<scitbx::vec3<double> >& xyz,
   const scitbx::af::const_ref<double>& boundary_layer_scaling_factors,
   const scitbx::af::const_ref<scitbx::vec3<double> >& h,
   const cctbx::xray::scattering_type_registry& registry) {
    scitbx::af::shared<std::complex<double> > amplitudes(h.size());
    std::string scattering_type;
    double f,f_bl,stol_sq,c,s;
    for (int i=0; i<h.size(); i++) {
      amplitudes[i] = std::complex<double>(0.0,0.0);
      stol_sq = 0.25*h[i]*h[i];
      f_bl = cctbx::eltbx::xray_scattering::wk1995("O",true).
        fetch().at_stol_sq(stol_sq);
      for (int j=0; j<scatterers.size(); j++) {
        scattering_type = scatterers[j];
        f = registry.gaussian(scattering_type)->at_stol_sq(stol_sq) +
          boundary_layer_scaling_factors[j]*f_bl;
        sincos(scitbx::constants::two_pi * (xyz[j]*h[i]),&s,&c);
        amplitudes[i] += std::complex<double>(f * c,f * s);
      }
    }
    return amplitudes;
  }

  /* ==========================================================================

     --------------------------------------------------------------------------
  */
  scitbx::af::shared<double> solvent_image
  (const scitbx::af::const_ref<double>& q) {
    scitbx::af::shared<double> s(q.size());
    double c = 1.0/scitbx::constants::two_pi;
    double h, stol_sq;
    cctbx::eltbx::xray_scattering::gaussian hoh =
      cctbx::eltbx::xray_scattering::wk1995("O",true).fetch();
    for (int i=0; i<s.size(); i++) {
      h = c * q[i];
      stol_sq = 0.25*h*h;
      s[i] = hoh.at_stol_sq(stol_sq);
      s[i] = s[i]*s[i];
    }

    return s;
  }

  /* ==========================================================================
     Finds the nearest neighbors

     --------------------------------------------------------------------------
  */
  scitbx::af::shared<int> nearest_neighbors
  (const scitbx::af::const_ref<scitbx::vec3<double> >& xyz,
   const scitbx::af::const_ref<double>& radius,
   const int& index, const double& probe_radius) {

    scitbx::af::shared<int> neighbors;
    scitbx::vec3<double> center = xyz[index];
    double base_distance = radius[index] + 2*probe_radius;
    double distance;
    for (int i=0; i<xyz.size(); i++) {
      distance = (xyz[i] - center).length();
      if (distance < (base_distance + radius[i])) {
        if (i != index) {
          neighbors.push_back(i);
        }
      }
    }

    return neighbors;
  }

  /* ==========================================================================
     Calculates the solvent accessible area
     Returns the fraction of the spherical area that is accessible where the
     radius of the sphere is the van der Waals radius
     --------------------------------------------------------------------------
  */
  scitbx::af::shared<int> solvent_accessible_area
  (const scitbx::af::const_ref<scitbx::vec3<double> >& xyz,
   const scitbx::af::const_ref<double>& radius,
   const scitbx::af::const_ref<int>& indices,
   const double& probe_radius=1.4,
   const int& n_points=1000) {

    scitbx::af::shared<int> areas(indices.size());
    scitbx::af::shared<scitbx::vec3<double> > base_sphere =
      bauer_spiral(n_points);
    scitbx::af::shared<int> neighbors;
    scitbx::vec3<double> sphere_point, current_center, neighbor_center;
    double current_radius, neighbor_radius, distance, n_accessible;
    bool point_is_accessible;
    for (int i=0; i<indices.size(); i++) {
      areas[i] = 0;
      n_accessible = 0.0;
      current_center = xyz[i];
      current_radius = radius[i] + probe_radius;
      neighbors = nearest_neighbors(xyz,radius,i,probe_radius);

      for (int j=0; j<base_sphere.size(); j++) {
        point_is_accessible = true;
        sphere_point = current_radius*base_sphere[j] + current_center;
        for (int k=0; k<neighbors.size(); k++) {
          neighbor_center = xyz[neighbors[k]];
          neighbor_radius = radius[neighbors[k]] + probe_radius;
          distance = (sphere_point - neighbor_center).length();

          if (distance < neighbor_radius) {
            point_is_accessible = false;
            break;
          }
        }
        if (point_is_accessible) {
          n_accessible += 1;
        }
      }

      areas[i] = n_accessible;
    }

    return areas;
  }

  /* ==========================================================================
     Applies a translation to a set of structure factors
     --------------------------------------------------------------------------
  */
  scitbx::af::shared<std::complex<double> > apply_translation
  (const scitbx::af::const_ref<std::complex<double> >& sf,
   const scitbx::af::const_ref<scitbx::vec3<double> >& h,
   const scitbx::vec3<double>& t) {

    assert(sf.size() == h.size());

    scitbx::af::shared<std::complex<double> > result(sf.size());
    double s,c;
    for (int i=0; i<sf.size(); i++) {
      sincos(scitbx::constants::two_pi * (h[i]*t),&s,&c);
      result[i] = sf[i] * std::complex<double>(c,s);
    }

    return result;
  }

  /* ==========================================================================
     c2 Class
     --------------------------------------------------------------------------
  */
  sastbx::fXS::c2::c2() {
    q_min = 0.0;
    q_max = 0.0;
    q_pixel_depth = 0;
    phi_pixel_radius = 0;
  }

  void sastbx::fXS::c2::set_image_composer(const sastbx::fXS::image_composer& a) {
    ic = a;
    q = ic.get_q();
    phi = ic.get_phi();
  }

  void sastbx::fXS::c2::set_q_limits(const double a, const double b) {
    q_min = a;
    q_max = b;
  }

  void sastbx::fXS::c2::set_pixel_sizes(const int a, const int b) {
    q_pixel_depth = a;
    phi_pixel_radius = b;
  }

  int sastbx::fXS::c2::get_q_index(const double& current_q) {
    int i_center = beam_center[1]*detector_size[0] + beam_center[0];
    int i_edge = i_center + (detector_size[0] - beam_center[0]) - 1;
    for (int i=i_center; i<i_edge; i++) {
      if (q[i] > current_q) {
        return i;
      }
    }
    return 0;
  }

  void sastbx::fXS::c2::initialize(bool verbose) {
    assert (q_min > 0.0);
    assert (q_max > q_min);
    assert (q_pixel_depth > 0);
    assert (phi_pixel_radius > 0);
    assert (q.size() > 0);
    assert (q.size() == phi.size());

    // make sure q_max is inside image
    beam_center = ic.get_beam_center();
    detector_size = ic.get_detector_size();
    int i_center = beam_center[1]*detector_size[0] + beam_center[0];
    int i_edge = i_center + (detector_size[0] - beam_center[0]) - 1;
    double max_q = q[i_edge];
    if (q_max > max_q) {
      q_max = max_q;
    }

    // determine number of rings and q limits for each ring
    int i_q_min = get_q_index(q_min);
    int i_q_max = get_q_index(q_max);
    scitbx::af::shared<double> ring_q_min;
    scitbx::af::shared<double> ring_q_max;
    double dq = q[i_center+1] - q[i_center];
    double q_i = q_min;
    while (q_i < q_max) {
      ring_q.push_back(q_i);
      ring_q_min.push_back(q_i - 0.5*q_pixel_depth*dq);
      ring_q_max.push_back(q_i + 0.5*q_pixel_depth*dq);
      q_i += (q_pixel_depth + 2)*dq;
    }

    // determine number of bins in each ring
    int current_i, n_bins;
    scitbx::af::shared<double> dphi( ring_q.size() );
    pixel_bins.resize( ring_q.size() );
    ring_c2.resize( ring_q.size() );
    for (int i=0; i<ring_q.size(); i++) {
      current_i = get_q_index(ring_q[i]);
      dphi[i] = phi[current_i + phi_pixel_radius*detector_size[0]] -
        (phi[current_i - phi_pixel_radius*detector_size[0]] -
         scitbx::constants::two_pi);
      n_bins = scitbx::math::ifloor(scitbx::constants::two_pi/dphi[i]);
      pixel_bins[i].resize(n_bins);
      ring_c2[i].resize(n_bins);
      dphi[i] = scitbx::constants::two_pi / n_bins;
    }

    // determine ring and bin position for each pixel
    int bin;
    for (int i=0; i<detector_size[0]*detector_size[1]; i++) {
      for (int j=0; j<ring_q.size(); j++) {
        if ( (q[i] >= ring_q_min[j]) && (q[i] < ring_q_max[j]) ) {
          bin = scitbx::math::ifloor(phi[i]/dphi[j]);
          pixel_bins[j][bin].push_back(i);
        }
      }
    }

    if (verbose) {
      for (int i=0; i<ring_q.size(); i++) {
        std::cout << i+1 << " " << ring_q[i] << " " << ring_c2[i].size() << "\n";
      }
    }
  }

  scitbx::af::shared<int> sastbx::fXS::c2::bin_mask() {
    scitbx::af::shared<int> bin_image(detector_size[0]*detector_size[1],0);
    int color;
    for (int i=0; i<pixel_bins.size(); i++) {
      for (int j=0; j<pixel_bins[i].size(); j++) {
        if (j%2 == 0) {
          color = 2;
        }
        else {
          color = 1;
        }
        for (int k=0; k<pixel_bins[i][j].size(); k++) {
          bin_image[pixel_bins[i][j][k]] = color;
        }
      }
    }
    return bin_image;
  }

  void sastbx::fXS::c2::process_intensities
  (const scitbx::af::const_ref<double>& intensities) {
    assert (intensities.size() == detector_size[0]*detector_size[1]);
    std::vector<double> mean_ring_intensities;
    double I_mean_sq;
    for (int i=0; i<pixel_bins.size(); i++) {            // loop over rings
      I_mean_sq = 0.0;

      // average intensities in bins
      mean_ring_intensities.clear();
      mean_ring_intensities.resize(pixel_bins[i].size(),0.0);
      for (int j=0; j<pixel_bins[i].size(); j++) {       // loop over bins
        for (int k=0; k<pixel_bins[i][j].size(); k++) {  // loop over pixels
          mean_ring_intensities[j] +=
            intensities[pixel_bins[i][j][k]];
        }
      }
      for (int j=0; j<pixel_bins[i].size(); j++) {
        mean_ring_intensities[j] =
          mean_ring_intensities[j]/pixel_bins[i][j].size();
      }

      // average intensity for ring
      for (int j=0; j<pixel_bins[i].size(); j++) {
        I_mean_sq += mean_ring_intensities[j];
      }
      I_mean_sq = I_mean_sq / mean_ring_intensities.size();
      I_mean_sq = I_mean_sq * I_mean_sq;

      // calculate autocorrelation function
      for (int shift=0; shift<mean_ring_intensities.size(); shift++) {
        ring_c2[i][shift] = 0.0;
        for (int j=0; j<mean_ring_intensities.size(); j++) {
          ring_c2[i][shift] +=
            ( mean_ring_intensities[j] *
              mean_ring_intensities[(j + shift)%mean_ring_intensities.size()] );
        }
        ring_c2[i][shift] = ring_c2[i][shift] / mean_ring_intensities.size() /
          I_mean_sq;
      }
    }
  }

  double sastbx::fXS::c2::get_ring_q(const int& i) {
    assert (i < ring_q.size());
    return ring_q[i];
  }

  int sastbx::fXS::c2::get_n_rings() {
    return ring_q.size();
  }

  scitbx::af::shared<double> sastbx::fXS::c2::get_c2(const int& i) {
    assert (i < ring_q.size());
    scitbx::af::shared<double> result(ring_c2[i].size());
    for (int j=0; j<result.size(); j++) {
      result[j] = ring_c2[i][j];
    }
    return result;
  }

  /* ==========================================================================
     old code for autocorrelation - bunch of functions
     --------------------------------------------------------------------------
  */
  scitbx::vec3<double> get_binning
  (const scitbx::vec2<int>& detector_size,
   const scitbx::vec2<int>& beam_center,
   const scitbx::af::const_ref<double>& q,
   const scitbx::af::const_ref<double>& phi,
   const double& ring_q) {

    // wedge is roughly 3 pixels deep, 4 pixels wide
    int q_pixel_depth = 3;
    int phi_pixel_radius = 2;
    int i_center = beam_center[1]*detector_size[0] + beam_center[0];
    double dq = q_pixel_depth*(q[i_center+1] - q[i_center]);
    int current_i;
    for (int i=i_center; i<(i_center + 0.5*detector_size[0]); i++) {
      if (q[i] > ring_q) {
        current_i = i;
        break;
      }
    }
    // does not work if pixel is off image
    double dphi = phi[current_i + phi_pixel_radius*detector_size[0]] -
      (phi[current_i - phi_pixel_radius*detector_size[0]] -
       scitbx::constants::two_pi);
    int n_bins = scitbx::math::ifloor(scitbx::constants::two_pi/dphi);
    dphi = scitbx::constants::two_pi / (n_bins);

    scitbx::vec3<double> binning;
    binning[0] = ring_q;
    binning[1] = dq;
    binning[2] = dphi;

    return binning;
  }

  scitbx::af::shared<int> prebin_image
  (const scitbx::vec3<double>& binning,
   const scitbx::af::const_ref<double>& q,
   const scitbx::af::const_ref<double>& phi) {

    int n_bins = scitbx::math::iceil(scitbx::constants::two_pi/binning[2]);
    scitbx::af::shared<int> bins(n_bins,0);
    double min_q = binning[0] - 0.5*binning[1];
    double max_q = binning[0] + 0.5*binning[1];
    int bin;
    for (int i=0; i<q.size(); i++) {
      if ( (q[i] >= min_q) && (q[i] < max_q) ) {
        bin = scitbx::math::ifloor(phi[i]/binning[2]);
        bins[bin] += 1;
      }
    }

    for (int i=0; i<n_bins; i++) {
      if (bins[i] == 0) {
        if (i == (n_bins - 1)) {
          bins.pop_back();
        }
        else {
          std::cout << "Warning: abnormal binning behavior\n"
                    << "bin " << i << " out of " << n_bins << " bins\n";
        }
      }
    }

    return bins;
  }

  scitbx::af::shared<int> bin_mask
  (const scitbx::vec3<double>& binning,
   const scitbx::af::const_ref<double>& q,
   const scitbx::af::const_ref<double>& phi) {

    scitbx::af::shared<int> bins(q.size(),1);
    double min_q = binning[0] - 0.5*binning[1];
    double max_q = binning[0] + 0.5*binning[1];
    int color = 1;
    for (int i=0; i<q.size(); i++) {
      if ( (q[i] >= min_q) && (q[i] < max_q) ) {
        if (scitbx::math::ifloor(phi[i]/binning[2])%2 == 0) {
          color = 3;
        }
        else {
          color = 2;
        }
        bins[i] = color;
      }
    }

    return bins;
  }

  scitbx::af::shared<double> bin_intensities
  (const scitbx::vec3<double>& binning,
   const scitbx::af::const_ref<double>& q,
   const scitbx::af::const_ref<double>& phi,
   const scitbx::af::const_ref<double>& intensities) {

    int n_bins = scitbx::math::iceil(scitbx::constants::two_pi/binning[2]);
    scitbx::af::shared<double> binned_intensities(n_bins,0.0);
    scitbx::af::shared<double> counts(n_bins,0.0);
    double min_q = binning[0] - 0.5*binning[1];
    double max_q = binning[0] + 0.5*binning[1];
    int bin;

    for (int i=0; i<q.size(); i++) {
      if ( (q[i] >= min_q) && (q[i] < max_q) ) {
        bin = scitbx::math::ifloor(phi[i]/binning[2]);
        binned_intensities[bin] += intensities[i];
        counts[bin] += 1.0;
      }
    }

    for (int i=0; i<n_bins; i++) {
      if (counts[i] == 0) {
        if (i == (n_bins - 1)) {
          binned_intensities.pop_back();
          counts.pop_back();
        }
        else {
          std::cout << "Warning: abnormal binning behavior\n"
                    << "bin " << i << " out of " << n_bins << " bins\n";
        }
      }
      else {
        binned_intensities[i] = binned_intensities[i]/counts[i];
      }
    }

    return binned_intensities;
  }

  /* ==========================================================================
     Equation 4 from Physical Review B 81, 174105 (2010) for one ring in one
     diffraction pattern
     --------------------------------------------------------------------------
  */
  scitbx::af::shared<double> autocorrelation
  (const scitbx::af::const_ref<double>& intensities) {
    scitbx::af::shared<double> ac(intensities.size());

    // square of average intensity
    double I_mean_sq = 0.0;
    for (int i=0; i<intensities.size(); i++) {
      I_mean_sq += intensities[i];
    }
    I_mean_sq = I_mean_sq / intensities.size();
    I_mean_sq = I_mean_sq * I_mean_sq;

    // scale autocorrelation by square of average intensity
    for (int shift=0; shift<intensities.size(); shift++) {
      ac[shift] = 0.0;
      for (int i=0; i<intensities.size(); i++) {
        ac[shift] += ( intensities[i] *
                       intensities[(i + shift)%intensities.size()] );
      }
      ac[shift] = ac[shift] / intensities.size() / I_mean_sq;
    }
    return ac;
  }

}
}
