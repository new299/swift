/*
    Swift (c) 2008 Genome Research Ltd.
    Authors: Nava Whiteford and Tom Skelly (new@sgenomics.org ts6@sanger.ac.uk)

    This file is part of Swift (http://swiftng.sourceforge.net).

    Swift is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Swift is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Swift.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SWIFTFFT_H_
#define SWIFTFFT_H_

#include <stdexcept>
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include "SwiftImage.h"
#include "SwiftImagePosition.h"
#include "Timetagger.h"
#include "stringify.h"

// Perform FFT on a SwiftImage. Cross-multiply with another SwiftFFT object
// and find the location of the maximum value. Uses FFTW (http://www.fftw.org)
// to do the heavy lifting. 
//
// This object is intended to be reusable to perform FFTs on multiple
// SwiftImages. That design avoids the need to repeatedly allocate and
// free a large (32MB for 2K*2K images) buffer, which can lead to
// fragmentation problems. It also amortises the overhead of using a
// FFTW_MEASURE plan, which takes some time to compute. (On that subject
// see the note in the constructor.)
//
// Unfortunately, reusability dictates that an image must first be loaded,
// then acted upon, in a series of method calls rather than all at once.
// This opens up the possibility that methods get called in the wrong order.
// The m_valid flag is used to indicate that a valid image has been loaded
// for subsequent processing by other methods.
//
// The SwiftImage contains pixel intensities stored as uint16. These
// need to be converted to doubles, stored in the 1-d buffer input to FFTW.
// 

using namespace std;

class SwiftFFT {

public:

  SwiftFFT (unsigned int image_width=0,
            unsigned int image_height=0,
            unsigned int fill_x=0,
            unsigned int fill_y=0,
            ostream &err_in=std::cerr)
          : m_image_width(image_width),
            m_image_height(image_height),
            m_fill_x(fill_x),
            m_fill_y(fill_y),
            m_image_valid(false),
            m_transform_valid(false),
            m_output_valid(false),
            m_offset (0,0) {

    // m_err << m_tt.str() << "Image size (x*y): " << m_image_width << " * " << m_image_height << endl;

    m_x_dim = pick_size (m_image_width+m_fill_x);
    m_y_dim = pick_size (m_image_height+m_fill_y);

    // m_err << m_tt.str() << "FFT size (x*y):   " << m_x_dim << " * " << m_y_dim << endl;

    // We'll allocate all the buffers we need right now. This saves a lot of freeing
    // and reallocating when the object is reused (as it is designed to be), and
    // gives us some exception safety as well, since all the frees are done in the
    // destructor.

    m_image_in   = (double*) fftw_malloc(sizeof(double) * m_x_dim * m_y_dim);
    m_magnitudes = (double*) fftw_malloc(sizeof(double) * m_x_dim * (m_y_dim/2+1));
    m_correl_out = (double*) fftw_malloc(sizeof(double) * m_x_dim * m_y_dim);
    m_fft_out    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_x_dim * (m_y_dim/2+1));
    m_cross      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_x_dim * (m_y_dim/2+1));

    // FFTW_MEASURE spends some time figuring out the best plan to use. FFTW_ESTIMATE
    // just guesses. Given our approach of sticking to simple batch sizes, it doesn't 
    // seem to matter much. For more complex batch sizes, it does -- see the comments
    // pertaining to the pick_size method.

    m_plan_forward = fftw_plan_dft_r2c_2d(m_x_dim, m_y_dim, m_image_in, m_fft_out, FFTW_MEASURE); 
    m_plan_reverse = fftw_plan_dft_c2r_2d(m_x_dim, m_y_dim, m_cross, m_correl_out, FFTW_MEASURE); 
    // m_err << m_tt.str() << "Forward and reverse FFT plans created" << endl;

  }
    
  ~SwiftFFT () {
        
    fftw_destroy_plan(m_plan_forward);
    fftw_destroy_plan(m_plan_reverse);
    fftw_free(m_image_in);
    fftw_free(m_magnitudes);
    fftw_free(m_correl_out);
    fftw_free(m_fft_out);
    fftw_free(m_cross);

  }

  int get_image_width ()  const {return m_image_width;}    
  int get_image_height () const {return m_image_height;}    

  fftw_complex *get_transformed_image () const {

    if ( ! m_transform_valid) {
      throw (std::domain_error("No transform done since image was loaded into SwiftFFT object"));
    }

    return m_fft_out;

  }

  SwiftImagePosition<> get_offset () const {

    if ( ! m_output_valid) {
      throw (std::domain_error("No cross-correlation done since image was loaded into SwiftFFT object"));
    }

    return m_offset;
  }

  // Load a SwiftImage to be transformed.
  //
  // Note that we're changing the index order here. A SwiftImage object stores the
  // image by rows -- the x coordinate varies the fastest, and so two consecutive
  // entries in the image vector are adjacent pixels in the same row. (That fact is
  // hidden within the operator() method.) FFTW expects the last index (y here) to
  // vary the fastest. And in the transformed image, which is complex, only y_max/2+1
  // points are stored -- i.e., the top half of the image. We could of course do the 
  // FFT "sideways", but that would be even more confusing than doing it this way.
  //
  // One further fiddle: According to the FFTW FAQ: "For human viewing of a spectrum,
  // it is often convenient to put the origin in frequency space at the center of the
  // output array, rather than in the zero-th element (the default in FFTW). If all
  // of the dimensions of your array are even, you can accomplish this by simply 
  // multiplying each element of the input array by (-1)^(i + j + ...), where
  // i, j, et cetera are the indices of the element. (This trick is a general property
  // of the DFT, and is not specific to FFTW.)" The even batchsize constraint is enforced
  // in pick_size. (NOTE: This feature has been removed!!!).
    
  bool load_SwiftImage (const SwiftImage<> &from_image, int at_x=0, int at_y=0) {
    int in_image_min_x = from_image.min_x();
    int in_image_max_x = from_image.max_x();
    int in_image_min_y = from_image.min_y();
    int in_image_max_y = from_image.max_y();

    int in_image_width  = from_image.max_x() - from_image.min_x();  
    int in_image_height = from_image.max_y() - from_image.min_y();  
    /*
    fprintf (stderr, "%sImage %d:%d (%d) %d:%d (%d)\n", m_tt.str(),
        in_image_min_x, in_image_max_x, in_image_width, 
        in_image_min_y, in_image_max_y, in_image_height); 
    */

    if (in_image_width  > m_image_width
        || in_image_height > m_image_height) {
      throw (std::domain_error("SwiftImage object input to SwiftFFT is too large."));
    }
    // Zero-fill unused areas.
    for (int ndx=0; ndx<m_x_dim*m_y_dim; ndx++) {
      m_image_in[ndx] = 0;
    }        

    // We want to map the origin of the input image to the point specified
    // by (at_x, at_y). Recall that the input image may have negative
    // coordinates. We will run the image-loading loop over the input image
    // coordinate range, and compute the correct porition for each input pixel.
    // We've already zeroed the entire array, so any array positions that don't
    // map to the input image will be zero.
    
    int to_x = at_x;

    for (int x=in_image_min_x; x<in_image_max_x; x++, to_x++) {

      int to_y = at_y;

      for (int y=in_image_min_y; y<in_image_max_y; y++, to_y++) {
        m_image_in[to_x*m_y_dim+to_y] = from_image(x,y);
      }

    }
    m_image_valid = true;
    m_transform_valid = false;
    m_output_valid = false;

    return true;

  }

  // Execute plan to create transform of loaded image.

  bool transform_image () {

    if ( ! m_image_valid) {
      throw (std::domain_error("No image has been loaded into SwiftFFT object"));
    }

    fftw_execute (m_plan_forward);
    m_transform_valid = 1;

    return true;

  }

  // Given another SwiftFFT object (the base image), perform FFT-driven cross-correlation
  // and update the offsets of this object to show offset relative to the base. This
  // just computes the cross-correlation (point-by-point multiply and sum) of the two images, 
  // but does it in the frequency domain because it's several orders of magnitude faster that
  // way than  using the direct approach. 

  bool cross_correlate (const SwiftFFT &other) {

    if ( ! m_transform_valid) {
      throw (std::domain_error("No transform done since image was loaded into SwiftFFT object"));
    }

    if (other.get_image_width() != m_image_width
        || other.get_image_height() != m_image_height) {
      throw (std::domain_error(string("Cannot cross-correlate differently-sized images: me: ") 
                                                       + stringify(m_image_width)   + string(",") + stringify(m_image_height) + 
                           string(" other: ") + stringify(other.get_image_width()) + string(",") + stringify(other.get_image_height())));
    }

    int y_out_dim = m_y_dim/2+1;             // y size of FFT output array

    fftw_complex *other_fft = other.get_transformed_image();

    // We're just multiplying two 1-d arrays point by point here, so a single
    // loop is enough.

    for (int ix=0; ix<m_x_dim*y_out_dim; ix++) {

      m_cross[ix][0] = m_fft_out[ix][0]*other_fft[ix][0]     // conjugate complex product
                     + m_fft_out[ix][1]*other_fft[ix][1];
      m_cross[ix][1] = m_fft_out[ix][1]*other_fft[ix][0]
                     - m_fft_out[ix][0]*other_fft[ix][1];

    }

    fftw_execute (m_plan_reverse);                            // reverse FFT of the result

    // Find the maximum. Its x,y offset is the image displacement. m_correl_out is
    // m_x_dim*m_y_dim real values. DC is at the corners. Also keep the second-best
    // peak, just for debugging purposes.

    double max_value = 0;
    int max_x = -1;
    int max_y = -1;

    double next_value = 0;
    int next_x = -1;
    int next_y = -1;

    for (int y=0; y<m_y_dim; y++) {
      for (int x=0; x<m_x_dim; x++) {

        int index = x*m_y_dim+y;

        if (m_correl_out[index] > max_value) { 
          
          next_value = max_value;                        // push old values down
          next_x = max_x;
          next_y = max_y;
          max_value = m_correl_out[index];
          max_x = x;
          max_y = y;
          
        } else if (m_correl_out[index] > next_value) {   // a new second-best?
          
          next_value = m_correl_out[index];
          next_x = x;
          next_y = y;
          
        }

      }
    }

    // Convert large positive offsets to small negative ones. Scale the result to account for 
    // zero-filling (including that which was done to get to a convenient FFT size, even if
    // fill_x and fill_y were zero).
    
    if (max_x > m_x_dim/2) {                // unwrap
//////      m_offset.x = ((max_x - m_x_dim) * m_image_width) / m_x_dim;
      m_offset.x = max_x - m_x_dim;
    } else {
//////      m_offset.x = (max_x * m_image_width) / m_x_dim;
      m_offset.x = max_x;
    }
    if (max_y > m_y_dim/2) {
//////      m_offset.y = ((max_y - m_y_dim) * m_image_height) / m_y_dim;
      m_offset.y = max_y - m_y_dim;
    } else {
//////      m_offset.y = (max_y * m_image_height) / m_y_dim;
      m_offset.y = max_y;
    }

    // Same for second-best. Offset values are local, not kept.
    
    int next_offset_x;
    int next_offset_y;
    
    if (next_x > m_x_dim/2) {                // unwrap
//////      next_offset_x = ((next_x - m_x_dim) * m_image_width) / m_x_dim;
      next_offset_x = next_x - m_x_dim;
    } else {
//////      next_offset_x = (next_x * m_image_width) / m_x_dim;
      next_offset_x = next_x;
    }
    if (next_y > m_y_dim/2) {
//////      next_offset_y = ((next_y - m_y_dim) * m_image_height) / m_y_dim;
      next_offset_y = next_y - m_y_dim;
    } else {
//////      next_offset_y = (next_y * m_image_height) / m_y_dim;
      next_offset_y = next_y;
    }

    m_output_valid = true;
/*
    cout << m_tt.str() << "Correl peak: " << max_value 
      << " at: " << m_offset.x << "  " << m_offset.y << endl;
    cout << m_tt.str() << "Correl 2nd:  " << next_value 
      << " at: " << next_offset_x << "  " << next_offset_y << endl;
*/
    return true;

  }
    
  // Return the magnitudes of the transformed image as a SwiftImage. The returned image
  // will be (slightly) larger than the original, since the FFT was padded out to a
  // convenient size for the FFT. It will also include the redundant, mirror-image
  // portion of the transform.
  //
  // Go read the note about index order in the header of load_SwiftImage. Then come back here.
  //
  // Besides having to flip the index order, we need to populate an entire image, but FFTW only
  // computes half (n/2+1, actually) that many points because the FFT of real data is conjugate
  // symmetrical. That is, for an M*N image, image[M-i, N-j] == image[i,j] conjugate.
  //
  // A new SwiftImage is created and returned to the caller, who is responsible for
  // destroying it.
  //
  // This method is intended for debugging. Actual Swift processing has no use for an image 
  // of the transformed input.

  SwiftImage<> *get_transformed_SwiftImage () {

    if ( ! m_transform_valid) {
      throw (std::domain_error("No transform done since image was loaded into SwiftFFT object"));
    }

    int y_out_dim = m_y_dim/2+1;             // y size of FFT output array

    // To scale the FFT output appropriately, we need to find the peak magnitude. We'll
    // need the amplitudes later, so save them. We won't use the DC component (near the
    // origin) as a maximum, since it tends to be larger than everything else. We'll 
    // truncate it later to max_uint16.

    double max_mag = 1;
    int max_x = -1;
    int max_y = -1;
    int avoid = 3;    // DC avoidance range

    for (int y=0; y<y_out_dim; y++) {
      for (int x=0; x<m_x_dim; x++) {

        int index = x*y_out_dim+y;
        double re = m_fft_out[index][0];
        double im = m_fft_out[index][1];
        double mag = sqrt(re*re + im*im);
        m_magnitudes[index] = mag;
        if (mag > max_mag 
            && abs(x-m_x_dim/2) > avoid
            && abs(y-m_y_dim/2) > avoid) {     // don't choose DC
          max_mag = mag;
            max_x = x;
            max_y = y;
        }

      }
    }

    cout << m_tt.str() << "Max magnitude: " << max_mag << " at: " << max_x << "  " << max_y << endl;

    double max_uint16 = 65535.0;
    double scaler = max_uint16 / max_mag;         // amount to scale uint16 FFT output by

    SwiftImage<> *si = new SwiftImage<> (m_x_dim, m_y_dim, 0);
    vector<uint16> &image = si->get_image();
    image.clear();               // vector was zero-filled

    for (int y=0; y<y_out_dim; y++) {             // e.g., y=0..32 for 64-point FFT
      for (int x=0; x<m_x_dim; x++) {
        int index = x*y_out_dim+y;
        double pixel = m_magnitudes[index] * scaler;
        if (pixel > max_uint16) {
          pixel = max_uint16;
        }
        image.push_back((uint16)pixel);
      }
    }

    // Second half: continue to load rows in increasing order, but load each row backwards.

    for (int y=y_out_dim; y<m_y_dim; y++) {        // y=33..63 for 64-point FFT
      for (int x=m_x_dim-1; x>= 0; x--) {        // load the row backwards
        int index = x*y_out_dim+m_y_dim-y;     // fetch from mirror-image row
        double pixel = m_magnitudes[index] * scaler;
        if (pixel > max_uint16) {
          pixel = max_uint16;
        }
        image.push_back((uint16)pixel);
      }
    }

    return si;

  }

  // Return the output of cross-correlation as a SwiftImage. The returned image
  // will be (slightly) larger than the original, since the FFT was padded out to a
  // convenient size for the FFT.
  //
  // A new SwiftImage is created and returned to the caller, who is responsible for
  // destroying it.
  //
  // This method is intended for debugging. Actual Swift processing has no use for an image 
  // of the cross-correlated input.

  SwiftImage<> *get_cross_correl_SwiftImage () {

    if ( ! m_output_valid) {
      throw (std::domain_error("No cross-correlation done since image was loaded into SwiftFFT object"));
    }

    // To scale the FFT output appropriately, we need to find the peak magnitude.

    int index = m_offset.x*m_y_dim + m_offset.y;
    double max_value = m_correl_out[index];
    double max_uint16 = 65535.0;
    double scaler = max_uint16 / max_value;         // amount to scale uint16 FFT output by

    SwiftImage<> *si = new SwiftImage<> (m_x_dim, m_y_dim, 0);
    vector<uint16> &image = si->get_image();
    image.clear();               // vector was zero-filled

    for (int y=0; y<m_y_dim; y++) {
      for (int x=0; x<m_x_dim; x++) {
        int index = x*m_y_dim+y;
        double pixel = m_correl_out[index] * scaler;
        if (pixel > max_uint16) {
          pixel = max_uint16;
        }
        image.push_back((uint16)pixel);
      }
    }

    return si;

  }

private:

  int          m_image_width;      // number of elements in a row
  int          m_image_height;     // number of elements in a column
  int          m_fill_x;           // x xero-fill amount
  int          m_fill_y;           // y zero-fill amount
  bool         m_image_valid;      // Has image been loaded?
  bool         m_transform_valid;  // Has transform been computed?
  bool         m_output_valid;     // Have offsets been computed?

  SwiftImagePosition<> m_offset;   // offset relative to reference image

  int          m_x_dim;            // FFT x batch size
  int          m_y_dim;            // FFT y batch size
  double       *m_image_in;        // pointer to input image data
  double       *m_magnitudes;      // pointer to saved magnitudes
  fftw_complex *m_fft_out;         // pointer to FFT output
  fftw_complex *m_cross;           // pointer to conjugate product of this with base image FFT
  double       *m_correl_out;      // pointer to reverse FFT of m_cross 

  fftw_plan    m_plan_forward;     // the cunning plan
  fftw_plan    m_plan_reverse;     // the cunning reverse plan

  Timetagger m_tt;

  static const int good_sizes[];   // table of good FFT batch sizes (see comment in .cpp)
  static const int num_sizes;      // size of above table
  int pick_size (int size);        // size-finding routine
  
};   // end of class definition

#endif /*SWIFTFFT_H_*/
