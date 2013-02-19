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

#ifndef SWIFTIMAGEANALYSIS_IMAGEANALYSIS_H
#define SWIFTIMAGEANALYSIS_IMAGEANALYSIS_H

#include "Cluster.h"
#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include "Timetagger.h"
#include "ChannelOffsets.h"
#include "SwiftImageCluster.h"
#include <math.h>
#include <string>
#include <algorithm>
#include <fstream>
#if defined(_OPENMP) 
#include <omp.h>
#endif


using namespace std;

template<class _prec=double,class _threshold_prec=uint16>
class ImageAnalysis {
public:
  
  typedef int base_type;            ///< type to use for bases

  static const int base_a       = 0; ///< These consts are the indexes which store the given base, this allows for accesses via base name, or iteration
  static const int base_c       = 1; ///< See base_a
  static const int base_g       = 2; ///< See base_a
  static const int base_t       = 3; ///< See base_a
  static const int base_invalid = 4; ///< See base_a

  static const int base_num     = 4; ///< Number of real bases (used to size image_filenames vector)


  ImageAnalysis(const vector<string> &image_filenames_a_in,       ///< List of filenames for A images, in cycle order
                const vector<string> &image_filenames_c_in,       ///< List of filenames for C images, in cycle order
                const vector<string> &image_filenames_g_in,       ///< List of filenames for G images, in cycle order
                const vector<string> &image_filenames_t_in,       ///< List of filenames for T images, in cycle order
                ostream &err_in=std::cerr                         ///< Write debug/errors here
                ) 
                : image_filenames(base_num),
                  images(base_num),
                  err(err_in) {

    // Add image filenames
    image_filenames[base_a] = image_filenames_a_in;
    image_filenames[base_c] = image_filenames_c_in;
    image_filenames[base_g] = image_filenames_g_in;
    image_filenames[base_t] = image_filenames_t_in;
  
    initialise();                
  }

  /// This constructor takes list a list of files, containing lists of files, which contain the images...
  ImageAnalysis(const string &image_filelist_filename_a,                            ///< File that contains a list of A image files in cycle order
                const string &image_filelist_filename_c,                            ///< File that contains a list of C image files in cycle order
                const string &image_filelist_filename_g,                            ///< File that contains a list of G image files in cycle order
                const string &image_filelist_filename_t,                            ///< File that contains a list of T image files in cycle order
                ostream &err_in=std::cerr                                           ///< Write debug/errors here
                ) 
                : image_filenames(base_num),
                  images(base_num),
                  err(err_in) {
  
    image_filenames[base_a] = read_image_list(image_filelist_filename_a);
    image_filenames[base_c] = read_image_list(image_filelist_filename_c);
    image_filenames[base_g] = read_image_list(image_filelist_filename_g);
    image_filenames[base_t] = read_image_list(image_filelist_filename_t);

    initialise();
  }
  
  void inline initialise();                 ///< Initialise this object, called by constructor

  void generate(vector<Cluster<_prec> > &clusters);              ///< Runs image analysis and returns a vector of clusters single pass analysis

public:
  
  int     params_threshold_window;            ///< Window size for thresholding, used to identify clusters
  _prec   params_threshold;                   ///< Thresholding parameter for image analysis, used to identift clusters
  bool    params_watershed;                   ///< Apply watershed segmentation?
  
  int     params_correlation_threshold_window;///< Window for image thresholding used in correlation (offset calculation)
  _prec   params_correlation_threshold;       ///< Threshold value used for offset calculation
  int     params_correlation_subimages;       ///< Number of subimages to calculate offsets for
  int     params_correlation_subsubimages;    ///< Number of subimages of subimages to use in calculation (the median of these is taken)
  int     params_correlation_cc_subimage_multiplier; ///< Multiply subimage by this when performing cross-channel offseting, allows crosschannel offset to be at a higher resolution.
  int     params_correlation_reference_cycle; ///< Reference cycle for offset calculation
  int     params_correlation_aggregate_cycle; ///< Aggregate cycle for offset calculation
  bool    params_correlation_median_channels; ///< Use median of different channels as offset?
  int     params_correlation_use_bases;       ///< Number of bases to use in channel offseting (default use all, -1)

  bool    params_background_subtraction_enabled; ///< Background subtraction in enabled?
  int     params_background_subtraction_window;///< Background subtraction window

  int     params_segment_cycles;              ///< Segment up to this many cycles

  bool    params_dump_images;                 ///< Dump images?

  bool    params_crop;                        ///< Crop input images (useful for removing flowcell wall)
  int     params_crop_start_x;                ///< Crop start x position
  int     params_crop_end_x;                  ///< Crop end x position
  int     params_crop_start_y;                ///< Crop start y position
  int     params_crop_end_y;                  ///< Crop end y position
 
  int     params_cluster_limit;               ///< maximum allowable number of clusters

  bool    params_remove_blended;              ///< Attempt to remove clusters still blended after deblending

  bool    params_calculate_noise;             ///< Attempt to calculate noise

  int     params_load_cycle;                   ///< Process this many cycles at a time

  ChannelOffsets<uint16>::correlation_type params_correlation_method;          ///< Image offset calculation method (not used)
  ChannelOffsets<_threshold_prec> *channel_offsets_thresholded;
  ChannelOffsets<uint16> *channel_offsets_standard;

  string  offsets_xml;                        ///< After a complete run this contains XML representing the offset maps applied to the images

private:

  void                         correlate_and_background_subtract();
  vector<SwiftImageCluster<_prec> > segmentation_and_registration(const vector<vector<SwiftImage<> > > &images);
  void                         build_clusters     (vector<Cluster<_prec> > &clusters, vector<SwiftImageCluster<_prec> > &image_clusters);
  void                         extend_clusters    (vector<Cluster<_prec> > &clusters, vector<SwiftImageCluster<_prec> > &image_clusters);
  void                         generate_initial   (vector<Cluster<_prec> > &clusters, vector<SwiftImageCluster<_prec> > &image_clusters);
  void                         generate_additional(vector<Cluster<_prec> > &clusters, vector<SwiftImageCluster<_prec> > &image_clusters);


  bool load_images(int load_first,bool grab_reference);                                             ///< Loads images into the images vector
  vector<string> read_image_list(string image_filelist_filename); ///< loads a file containing filelists and returns it as a string
  
  int total_cycles;                           ///< total number of cycles, populated by load_images

  vector<vector<string> > image_filenames;    ///< vector, of filename vectors Indexed using base_a/t/g/c
  vector<vector<SwiftImage<uint16> > > images;///< this vector holds the actual image data, it is populated by load_images

  vector<SwiftImage<uint16> >          reference_images;
  
  ostream &err;                               ///< Error output will be writen here, set to cerr in constructor default
  Timetagger m_tt;                            ///< Timetag-generating object
  
};

// This needs to be included here because this is a template class
#include "ImageAnalysis.cpp"

#endif
