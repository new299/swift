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

#ifndef SWIFT_MOCKIMAGEANALYSIS_H
#define SWIFT_MOCKIMAGEANALYSIS_H

#include "Cluster.h"
#include <iostream>
#include <vector>

using namespace std;

template<class _prec=double>
class MockImageAnalysis {
public:

  static const int base_a       = 0; ///< These consts are the indexes which store the given base, this allows for accesses via base name, or iteration
  static const int base_c       = 1; ///< See base_a
  static const int base_g       = 2; ///< See base_a
  static const int base_t       = 3; ///< See base_a
  static const int base_invalid = 4; ///< See base_a

  static const int base_num     = 4; ///< Number of real bases (used to size image_filenames vector)

  MockImageAnalysis(const vector<string> &image_filenames_a_in,       ///< List of filenames for A images, in cycle order
                const vector<string> &image_filenames_t_in,       ///< List of filenames for T images, in cycle order
                const vector<string> &image_filenames_g_in,       ///< List of filenames for G images, in cycle order
                const vector<string> &image_filenames_c_in,       ///< List of filenames for C images, in cycle order
                _prec fwhm_in = default_cluster_fwhm,                               ///< FWHM input parameter
                _prec distance_in = default_cluster_distance,                       ///< Min Cluster distance 
                _prec threshold_in = default_image_threshold,                       ///< Image thresholding
                _prec image_max_stage_offset_in = default_image_max_stage_offset,   ///< Image maximum offset
                _prec image_max_colour_offset_in = default_image_max_colour_offset, ///< Image colour offset?
                ostream &err_in=std::cerr                                           ///< Write debug/errors here
                ) 
                : cluster_fwhm(fwhm_in),
                  cluster_distance(distance_in),
		  image_threshold(threshold_in),
                  image_max_stage_offset(image_max_stage_offset_in),
                  image_max_colour_offset(image_max_colour_offset_in),
		  image_filenames(base_num),
                  err(err_in) {

    // Add image filenames
    image_filenames[base_a] = image_filenames_a_in;
    image_filenames[base_t] = image_filenames_t_in;
    image_filenames[base_g] = image_filenames_g_in;
    image_filenames[base_c] = image_filenames_c_in;
  
    initialise();                
  }

  /// This constructor takes list a list of files, containing lists of files, which contain the images...
  MockImageAnalysis(const string &image_filelist_filename_a,                            ///< File that contains a list of A image files in cycle order
                const string &image_filelist_filename_t,                            ///< File that contains a list of T image files in cycle order
                const string &image_filelist_filename_g,                            ///< File that contains a list of G image files in cycle order
                const string &image_filelist_filename_c,                            ///< File that contains a list of C image files in cycle order
                _prec fwhm_in = default_cluster_fwhm,                               ///< FWHM input parameter
                _prec distance_in = default_cluster_distance,                       ///< Min Cluster distance 
                _prec threshold_in = default_image_threshold,                       ///< Image thresholding
                _prec image_max_stage_offset_in = default_image_max_stage_offset,   ///< Image maximum offset
                _prec image_max_colour_offset_in = default_image_max_colour_offset, ///< Image colour offset?
                ostream &err_in=std::cerr                                           ///< Write debug/errors here
                ) 
                : cluster_fwhm(fwhm_in),
                  cluster_distance(distance_in),
		  image_threshold(threshold_in),
                  image_max_stage_offset(image_max_stage_offset_in),
                  image_max_colour_offset(image_max_colour_offset_in),
		  image_filenames(base_num),
                  err(err_in) {
  
    image_filenames[base_a] = read_image_list(image_filelist_filename_a);
    image_filenames[base_t] = read_image_list(image_filelist_filename_t);
    image_filenames[base_g] = read_image_list(image_filelist_filename_g);
    image_filenames[base_c] = read_image_list(image_filelist_filename_c);

    initialise();
  }

  void inline initialise();                 ///< Initialise this object, called by constructor

  bool run_image_analysis(int pass=0); ///< Runs the image analysis

  inline const vector<Cluster<_prec> > &generate();              ///< Runs image analysis and returns a vector of clusters single pass analysis
  inline const vector<Cluster<_prec> > &generate_1pass();        ///< Runs image analysis and returns a vector of clusters single pass analysis
  inline const vector<Cluster<_prec> > &generate_2pass();        ///< Runs image analysis and returns a vector of clusters two pass analysis

  void average_offsets();                                   ///< Currently calls Klaus's python script to average image offsets (really messy external call)
  void save_alignment_transforms(string filename);          ///< Saves all alignment transforms in format Klaus's code can read
  void load_default_alignment_transforms(string filename);  ///< Loads alignment transforms, as produced by Klaus's code

  const vector<Cluster<_prec> > &get_clusters() { ///< Returns clusters that have been constructed during image analysis
    return clusters;
  }

private:
  
  bool load_images();                                             ///< Loads images into the images vector
  vector<string> read_image_list(string image_filelist_filename); ///< loads a file containing filelists and returns it as a string
  
  int total_cycles;                           ///< total number of cycles, populated by load_images
  _prec cluster_fwhm;                         ///< Cluster fwhm, analysis input parameter
  _prec cluster_distance;                     ///< Cluster min distance, analysis input parameter
  _prec image_threshold;                      ///< Image threshold, analysis input parameter
  _prec image_max_stage_offset;               ///< Image maximum offset
  _prec image_max_colour_offset;              ///< Image colour offset

  vector<vector<string> > image_filenames;    ///< vector, of filename vectors Indexed using base_a/t/g/c
  vector<Cluster<_prec> > clusters;                ///< Stores resulting clusters

  ostream &err;                                                ///< Error output will be writen here, set to cerr in constructor default
  
  static const _prec default_cluster_fwhm            = 2.7;    ///< default value for cluster? fwhm. 
  static const _prec default_cluster_distance        = 1.5;    ///< default value for cluster min distance. 
  static const _prec default_image_threshold         = 4.0;    ///< default value for image thresholding.
  static const _prec default_image_max_stage_offset  = 500.0;  ///< default value for image max stage offset
  static const _prec default_image_max_colour_offset = 5.0;    ///< default value for image colour offset
};

// This needs to be included here because this is a template class
#include "MockImageAnalysis.cpp"

#endif
