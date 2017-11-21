#include <sastbx/fXS/detector_tools.h>

namespace sastbx {
namespace fXS {

  /* ==========================================================================
     detector_geometry Class

     Parameters describing detector (units must be consistent)
     Parameters follow the conventions established for the CXI format
     (http://www.cxidb.org/cxi.html)
     --------------------------------------------------------------------------
  */

  sastbx::fXS::detector_geometry::detector_geometry() {}

  // corner position is (x, y, z) in units of length
  void sastbx::fXS::detector_geometry::set_corner_position
  (const scitbx::vec3<double>& cp) {
    corner_position = cp;
  }

  scitbx::vec3<double> sastbx::fXS::detector_geometry::get_corner_position() {
    return corner_position;
  }

  // detector size is (n_rows, n_columns) in pixels
  void sastbx::fXS::detector_geometry::set_detector_size
  (const scitbx::vec2<int>& ds) {
    detector_size = ds;
  }

  scitbx::vec2<int> sastbx::fXS::detector_geometry::get_detector_size() {
    return detector_size;
  }

  // pixel size is (size along x, size along y) in units of length
  void sastbx::fXS::detector_geometry::set_pixel_size
  (const scitbx::vec2<double>& ps) {
    pixel_size = ps;
  }

  scitbx::vec2<double> sastbx::fXS::detector_geometry::get_pixel_size() {
    return pixel_size;
  }

  /* ==========================================================================
     image_base Class

     Takes experimental parameters and converts them into reciprocal space
     parameters (h and q), and phi
     --------------------------------------------------------------------------
  */

  sastbx::fXS::image_base::image_base() {}

  void sastbx::fXS::image_base::set_ewald_sphere(const ewald_sphere& esa) {
    es = esa;
  }

  ewald_sphere sastbx::fXS::image_base::get_ewald_sphere() {
    return es;
  }

  void sastbx::fXS::image_base::set_detector_geometry
  (const detector_geometry& dga) {
    dg = dga;
  }

  detector_geometry sastbx::fXS::image_base::get_detector_geometry() {
    return dg;
  }

  // reset q, phi, and h
  void sastbx::fXS::image_base::reset() {
    q.clear();
    phi.clear();
    h.clear();
  }

  // q values at center of each pixel
  scitbx::af::shared<double> sastbx::fXS::image_base::get_center_q() {
    if (q.size() == 0) {
      scitbx::vec3<double> corner_position = dg.get_corner_position();
      scitbx::vec2<int> detector_size = dg.get_detector_size();
      scitbx::vec2<double> pixel_size = dg.get_pixel_size();
      scitbx::vec3<double> current_h;
      scitbx::vec2<double> current_xy;
      double current_phi;
      q.reserve( detector_size[0] * detector_size[1] );
      phi.reserve( detector_size[0] * detector_size[1] );
      es.set_distance(corner_position[2]);
      for (int row=0; row<detector_size[0]; row++) {
        for (int col=0; col<detector_size[1]; col++) {
          current_xy[0] = corner_position[0] - (col + 0.5)*pixel_size[0];
          current_xy[1] = corner_position[1] - (row + 0.5)*pixel_size[1];
          current_h = es.get_h(current_xy);
          q.push_back(es.h_to_q(current_h));
          current_phi = std::atan(current_xy[1]/current_xy[0]);
          if (current_phi < 0) {
            current_phi += scitbx::constants::pi;
          }
          if (current_xy[1] < 0.0) {
            current_phi += scitbx::constants::pi;
          }
          phi.push_back(current_phi);
        }
      }
    }
    return q;
  }

  // phi values at center of each pixel
  scitbx::af::shared<double> sastbx::fXS::image_base::get_center_phi() {
    if (phi.size() == 0) {
      get_center_q();
    }
    return phi;
  }

  // h values at the corners of each pixel
  scitbx::af::shared<scitbx::vec3<double> >
  sastbx::fXS::image_base::get_corner_h() {
    if (h.size() == 0) {
      scitbx::vec3<double> corner_position = dg.get_corner_position();
      scitbx::vec2<int> detector_size = dg.get_detector_size();
      scitbx::vec2<double> pixel_size = dg.get_pixel_size();
      int n_row = detector_size[0] + 1;
      int n_col = detector_size[1] + 1;
      scitbx::vec3<double> current_h;
      scitbx::vec2<double> current_xy;
      int current_index;
      h.reserve( n_row * n_col );
      es.set_distance(corner_position[2]);
      for (int row=0; row<n_row; row++) {
        for (int col=0; col<n_col; col++) {
          current_index = row*detector_size[1] + col;
          current_xy[0] = corner_position[0] - col*pixel_size[0];
          current_xy[1] = corner_position[1] - row*pixel_size[1];
          current_h = es.get_h(current_xy);
          h.push_back(current_h);
        }
      }
    }
    return h;
  }

  // integrate intensity over area of each pixel
  scitbx::af::shared<double>
  sastbx::fXS::image_base::integrate
  (const scitbx::af::const_ref<double>& intensities) {

    assert(intensities.size() == h.size());
    scitbx::vec2<int> detector_size = dg.get_detector_size();
    scitbx::vec2<double> pixel_size = dg.get_pixel_size();
    int dsx = detector_size[0] + 1;
    scitbx::af::shared<double> pixels(detector_size[0]*detector_size[1]);
    double I00,I01,I10,I11, pixel_area;
    pixel_area = pixel_size[0]*pixel_size[1];
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

  /* ==========================================================================
     c2_tiles Class

     Handles the autocorrelation calculation for detectors composed of multiple
     tiles
     --------------------------------------------------------------------------
  */

  sastbx::fXS::c2_tiles::c2_tiles() {
    upper_left_xy[0] = -1.0e12;
    upper_left_xy[1] = -1.0e12;
    lower_right_xy[0] = 1.0e12;
    lower_right_xy[1] = 1.0e12;
  }

  void sastbx::fXS::c2_tiles::set_ewald_sphere(const ewald_sphere& esa) {
    ib.set_ewald_sphere(esa);
  }

  void sastbx::fXS::c2_tiles::add_geometry(const detector_geometry& dg) {
    tiles.push_back(dg);
    ib.set_detector_geometry(dg);
    ib.reset();
    q.push_back(ib.get_center_q().deep_copy());
    phi.push_back(ib.get_center_phi().deep_copy());
    scitbx::vec3<double> cp = ib.get_detector_geometry().get_corner_position();
    scitbx::vec2<double> ps = ib.get_detector_geometry().get_pixel_size();
    scitbx::vec2<int> ds = ib.get_detector_geometry().get_detector_size();
    double tmp;
    for (int i=0; i<2; i++) {
      if (cp[i] > upper_left_xy[i]) {
        upper_left_xy[i] = cp[i];
      }
    }
    tmp = cp[0] - ps[1]*ds[1];
    if ( tmp  < lower_right_xy[0] ) {
      lower_right_xy[0] = tmp;
    }
    tmp = cp[1] - ps[0]*ds[0];
    if ( tmp  < lower_right_xy[1] ) {
      lower_right_xy[1] = tmp;
    }
  }

  void sastbx::fXS::c2_tiles::add_intensities
  (const scitbx::af::const_ref<double>& i_data) {
    intensities.push_back(i_data);
  }

  void sastbx::fXS::c2_tiles::set_q_limits(const double& a, const double& b) {
    q_min = a;
    q_max = b;
  }

  void sastbx::fXS::c2_tiles::set_ring_pixel_sizes(const int& a, const int& b) {
    q_pixel_depth = a;
    phi_pixel_radius = b;
  }

  int sastbx::fXS::c2_tiles::get_i_from_q0(const double& q_i) {
    int current_index;
    for (int i=0; i<q0.size(); i++) {
      if (q_i > q0[i]) {
        current_index = i;
        break;
      }
    }
    return current_index;
  }

  void sastbx::fXS::c2_tiles::initialize() {

    // check parameters
    scitbx::math::basic_statistics<double> bs;
    assert (tiles.size() > 0);
    assert (q_max > q_min);
    assert (q_pixel_depth > 0);
    assert (phi_pixel_radius > 0);
    double min_q = 1.0e12;
    double max_q = -1.0e12;
    // assume all tiles have same pixel sizes
    scitbx::vec2<double> pixel_size = tiles[0].get_pixel_size();
    for (int i=0; i<tiles.size(); i++) {
      assert (q[i].size() > 0);
      bs = scitbx::math::basic_statistics<double>(q[i].const_ref());
      if (bs.min < min_q) {
        min_q = bs.min;
      }
      if (bs.max > max_q) {
        max_q = bs.max;
      }
    }

    // get q values along (x,0,z) to calculate dq per pixel
    detector_geometry tmp_dg;
    scitbx::vec3<double> tmp_cp = tiles[0].get_corner_position();
    tmp_cp[0] = 2*upper_left_xy[0];
    tmp_cp[1] = 0.0;
    tmp_dg.set_corner_position(tmp_cp);
    scitbx::vec2<int> tmp_ds;
    tmp_ds[0] = phi_pixel_radius;
    tmp_ds[1] = 2*int((upper_left_xy[0] - lower_right_xy[1])/pixel_size[0]) + 1;
    tmp_dg.set_detector_size(tmp_ds);
    tmp_dg.set_pixel_size(pixel_size);
    ib.set_detector_geometry(tmp_dg);
    ib.reset();
    q0 = ib.get_center_q();
    phi0 = ib.get_center_phi();
    double dq = std::fabs(q0[0] - q0[1]);
    if (q_min < min_q) {
      q_min = min_q + dq*q_pixel_depth;
    }
    if ( (q_max > max_q) || (q_max < q_min) ) {
      q_max = max_q - dq*q_pixel_depth;
    }

    // determine number of rings and q limits for each ring
    scitbx::af::shared<double> ring_q_min;
    scitbx::af::shared<double> ring_q_max;
    double q_i = q_min;
    while (q_i < q_max) {
      ring_q.push_back(q_i);
      ring_q_min.push_back(q_i - 0.5*dq*q_pixel_depth);
      ring_q_max.push_back(q_i + 0.5*dq*q_pixel_depth);
      q_i += (q_pixel_depth + 2)*dq;
    }

    // determine number of bins in each ring
    int n_bins, i_q_ring;
    scitbx::af::shared<double> dphi( ring_q.size() );
    pixel_bins.resize( ring_q.size() );
    mean_ring_intensities.resize( ring_q.size() );
    ring_intensities.resize( ring_q.size() );
    ring_c2.resize( ring_q.size() );
    for (int i=0; i<ring_q.size(); i++) {
      i_q_ring = get_i_from_q0(ring_q[i]);
      dphi[i] = ( scitbx::constants::two_pi -
                  phi0[(phi_pixel_radius - 1)*tmp_ds[1] + i_q_ring] );
      n_bins = scitbx::math::ifloor(scitbx::constants::two_pi/dphi[i]);
      dphi[i] = scitbx::constants::two_pi/n_bins;
      pixel_bins[i].resize(n_bins);
      for (int j=0; j<n_bins; j++) {
        pixel_bins[i][j].resize(tiles.size());
      }
      ring_c2[i].resize(n_bins);
      dphi[i] = scitbx::constants::two_pi / n_bins;
    }

    // determine ring and bin position for each pixel
    int bin;
    for (int tile=0; tile<tiles.size(); tile++) {
      tmp_ds = tiles[tile].get_detector_size();
      for (int i=0; i<tmp_ds[0]*tmp_ds[1]; i++) {
        for (int j=0; j<ring_q.size(); j++) {
          if ((q[tile][i] >= ring_q_min[j]) && (q[tile][i] < ring_q_max[j])) {
            bin = scitbx::math::ifloor(phi[tile][i]/dphi[j]);
            pixel_bins[j][bin][tile].push_back(i);
          }
        }
      }
    }

    // check for bins that are at least half full
    //int cutoff = q_pixel_depth * phi_pixel_radius / 2;
    int cutoff = 0;
    use_bin.resize( pixel_bins.size() );
    for (int i=0; i<pixel_bins.size(); i++) {
      use_bin[i].resize( pixel_bins[i].size() );
      for (int j=0; j<pixel_bins[i].size(); j++) {
        use_bin[i][j].resize( pixel_bins[i][j].size() );
        for (int k=0; k<pixel_bins[i][j].size(); k++) {
          if (pixel_bins[i][j][k].size() > cutoff) {
            use_bin[i][j][k] = true;
          }
          else {
            use_bin[i][j][k] = false;
          }
        }
      }
    }

  }

  scitbx::af::shared<int> sastbx::fXS::c2_tiles::bin_mask(const int& tile) {
    scitbx::vec2<int> ds = tiles[tile].get_detector_size();
    scitbx::af::shared<int> bin_image(ds[0]*ds[1],0);
    int color;
    for (int i=0; i<pixel_bins.size(); i++) {
      for (int j=0; j<pixel_bins[i].size(); j++) {
        if (j%2 == 0) {
          color = 2;
        }
        else {
          color = 1;
        }
        if (use_bin[i][j][tile]) {
          for (int k=0; k<pixel_bins[i][j][tile].size(); k++) {
            bin_image[pixel_bins[i][j][tile][k]] = color;
          }
        }
      }
    }
    return bin_image;
  }
  void sastbx::fXS::c2_tiles::process_intensities() {
    assert (intensities.size() == tiles.size());
    for (int i=0; i<pixel_bins.size(); i++) {            // loop over rings
      process_ring(i);
    }
  }

  void sastbx::fXS::c2_tiles::process_ring(const int& i) {
    assert (intensities.size() == tiles.size());
    std::vector<double> current_ring_intensities(pixel_bins[i].size(),0.0);
    double I_mean_sq = 0.0;
    int bin_count = 0;
    int ring_pixel_count = 0;
    int bin_pixel_count;

    // average intensities in bins
    for (int j=0; j<pixel_bins[i].size(); j++) {       // loop over bins
      bin_pixel_count = 0;
      for (int tile=0; tile<pixel_bins[i][j].size(); tile++) { // loop over tiles
        if (use_bin[i][j][tile]) {
          for (int k=0; k<pixel_bins[i][j][tile].size(); k++) {// loop over pixels
            current_ring_intensities[j] +=
              intensities[tile][pixel_bins[i][j][tile][k]];
            I_mean_sq += intensities[tile][pixel_bins[i][j][tile][k]];
          }
          bin_pixel_count = pixel_bins[i][j][tile].size();
          ring_pixel_count += bin_pixel_count;
          bin_count += 1;
        }
      }
      if (bin_pixel_count != 0) {
        current_ring_intensities[j] = current_ring_intensities[j]/bin_pixel_count;
      }
    }

    // average intensity for ring
    I_mean_sq = I_mean_sq / ring_pixel_count;
    mean_ring_intensities[i] = I_mean_sq;
    ring_intensities[i] = current_ring_intensities;
    I_mean_sq = I_mean_sq * I_mean_sq;

    // calculate autocorrelation function
    for (int shift=0; shift<current_ring_intensities.size(); shift++) {
      ring_c2[i][shift] = 0.0;
      for (int j=0; j<current_ring_intensities.size(); j++) {
        ring_c2[i][shift] +=
          ( current_ring_intensities[j] *
            current_ring_intensities[(j + shift)%current_ring_intensities.size()] );
      }
      ring_c2[i][shift] = ring_c2[i][shift] / bin_count / I_mean_sq;
    }
  }

  void sastbx::fXS::c2_tiles::reset_intensities() {
    intensities.clear();
  }

  double sastbx::fXS::c2_tiles::get_ring_q(const int& i) {
    assert (i < ring_q.size());
    return ring_q[i];
  }

  scitbx::af::shared<double> sastbx::fXS::c2_tiles::get_mean_ring_intensities() {
    return mean_ring_intensities;
  }

  scitbx::af::shared<double> sastbx::fXS::c2_tiles::get_ring_intensities
  (const int& i) {
    assert (i <ring_q.size());
    scitbx::af::shared<double> result(ring_intensities[i].size());
    for (int j=0; j<ring_intensities[i].size(); j++) {
      result[j] = ring_intensities[i][j];
    }
    return result;
  }

  int sastbx::fXS::c2_tiles::get_n_rings() {
    return ring_q.size();
  }

  scitbx::af::shared<double> sastbx::fXS::c2_tiles::get_c2(const int& i) {
    assert (i < ring_q.size());
    scitbx::af::shared<double> result(ring_c2[i].size());
    for (int j=0; j<result.size(); j++) {
      result[j] = ring_c2[i][j];
    }
    return result;
  }

  scitbx::af::shared<int> sastbx::fXS::c2_tiles::get_pixel_indices
  (const int& i, const int& tile) {
    assert (i < ring_q.size());
    scitbx::vec2<int> tmp_ds = tiles[tile].get_detector_size();
    int row, col;
    std::set<int> pixel_indices;
    for (int j=0; j<pixel_bins[i].size(); j++) {
      if (use_bin[i][j][tile]) {
        for (int k=0; k<pixel_bins[i][j][tile].size(); k++) {
          row = pixel_bins[i][j][tile][k]/tmp_ds[1];
          col = pixel_bins[i][j][tile][k]%tmp_ds[1];
          pixel_indices.insert(row*(tmp_ds[1] + 1) + col);
          pixel_indices.insert(row*(tmp_ds[1] + 1) + col + 1);
          pixel_indices.insert((row + 1)*(tmp_ds[1] + 1) + col);
          pixel_indices.insert((row + 1)*(tmp_ds[1] + 1) + col + 1);
        }
      }
    }
    scitbx::af::shared<int> result;
    for (std::set<int>::iterator iter = pixel_indices.begin();
         iter != pixel_indices.end(); ++iter) {
      result.push_back(*iter);
    }
    return result;
  }

}
}
