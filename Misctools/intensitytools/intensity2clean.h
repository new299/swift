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

#ifndef INTENSITYTOOLS_INTENSITY2FLOWCELLREMOVE_H
#define INTENSITYTOOLS_INTENSITY2FLOWCELLREMOVE_H

#include "intensity2xybins.h"
#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;

/// This function detects a hump low to high to low transition in a series of bins
/// the transition threshold is the average*threshold
/// start < 
int find_hump(vector<int> &bins, ///< A vector containing the number of items in each bin
              double threshold,  ///< Scaler for average to create threshold for high/low transitions
              int start,         ///< Start processing from this bin
              int end,           ///< Stop at this bin
              int direction,     ///< Increment in this direction i.e. 1 = 0->1, -1 = 1->0
              ostream &out       ///< Write errors/debugging info to this stream
              ) {
  
  // find max
  double average_in_bins=0;
  for(vector<int>::iterator i = bins.begin();i != bins.end();i++) {
    average_in_bins += (*i);
  }
  average_in_bins = average_in_bins/bins.size();
  out << "average: " << average_in_bins << endl; 

  bool low=false;
  bool high=false;
  bool first=true;

  bool high_low_trans=false;
  bool low_high_trans=false;

  // Walk across bins, detect transition
  int hump_end=start;
  for(int n=start;(n != end) && (hump_end == 0);n=n+direction) {
    if(bins[n] < (average_in_bins*threshold)) {
      low=true;
      if(!first) {
        if(high) {high_low_trans = true;}
        high=false;
      }
      first=false;
    }
    else {
      high=true;
      if(!first) {
        if(low) {low_high_trans = true;}
        low=false;
      } 
      first=false;
   }
   if(low_high_trans && high_low_trans) {
     hump_end=n;
   }
  }

  return hump_end;
}


/// This function is designed to detect a flowcell wall in an intensity file and remove it
/// It does this by detecting humps in the distribution of clusters in the X direction.
/// It takes a set of input files (intensity,noise and positions) and outputs a set of
/// 'good' files and 'bad' files. 'Good' files contain intensities without flowcell walls
/// 'bad' files contain the clusters that were removed. The noise and position files are
/// not used in processing, they are simply kept in sync with the intensity files.
void intensity2flowcellremove(char *input_intensityfilename,             ///< Input intensity filename
                              char *input_noisefilename,                 ///< Input noise filename
                              char *input_positionsfilename,             ///< Input positions filename
                              char *output_good_intensityfilename,       ///< Output intensities with flowcell removed filename
                              char *output_good_noisefilename,           ///< Output noise with flowcell removed filename
                              char *output_good_positionsfilename,       ///< Output positions with flowcell removed filename
                              char *output_bad_intensityfilename,        ///< Output intensities that were removed filename
                              char *output_bad_noisefilename,            ///< Output noise that was removed filename
                              char *output_bad_positionsfilename,        ///< Output positions that were removed filename
                              int bin_count,                             ///< The number of bins to use in detection
                              double threshold,                          ///< The threshold scaler to use in detectino
                              ostream &out                               ///< The stream to write errors/debugging info to
                             ) {                            
  
  
  // Calculate bins
  vector<int> xbins;
  vector<int> ybins;
  double bin_size=0;
  intensity2xybins(input_intensityfilename,bin_count,bin_size,xbins,ybins);

  // find hump in xbins, a hump means the bins rise above 25% of the maximum value and then fall below 25% again within the first 25% of bins.
  
  int left_hump_end  = find_hump(xbins,threshold,0,xbins.size()/4,1,out);
  int right_hump_end = find_hump(xbins,threshold,xbins.size(),xbins.size()-(xbins.size()/4),-1,out);
 
  // out << "Left Hump End Bin: " << left_hump_end << endl;
  // out << "Left Hump End Pos: " << left_hump_end*bin_size << endl;
  
  // out << "Right Hump End Bin: " << right_hump_end << endl;
  // out << "Right Hump End Pos: " << right_hump_end*bin_size << endl;

  // Push these out a little to avoid numerical errors removing some reads.
  if(left_hump_end == 0) left_hump_end = -1;
  if(right_hump_end == bin_count) right_hump_end = bin_count+1;
  
  out << "Using Left Hump End Bin: " << left_hump_end << endl;
  out << "Using Left Hump End Pos: " << left_hump_end*bin_size << endl;
  
  out << "Using Right Hump End Bin: " << right_hump_end << endl;
  out << "Using Right Hump End Pos: " << right_hump_end*bin_size << endl;
  
  
  // Now we had identified the hump end position delete all reads before the hump

  // Intensity filename reader
  IntensityReader input_intensity_file (input_intensityfilename);
  IntensityReader input_noise_file     (input_noisefilename);
  ifstream input_positions_file(input_positionsfilename);
  
  IntensityWriter output_good_intensity_file(output_good_intensityfilename);
  IntensityWriter output_good_noise_file    (output_good_noisefilename);
  ofstream        output_good_positions_file(output_good_positionsfilename);
  
  IntensityWriter output_bad_intensity_file;
  IntensityWriter output_bad_noise_file;
  ofstream        output_bad_positions_file;

  // Filter clusters based on the humps we identified.
  int discarded=0;
  for(;!input_intensity_file.eof();) {
   
    // An intensitysequence is a single line of the intensity filename
    IntensitySequence is = input_intensity_file.get_next();
    IntensitySequence in = input_noise_file.get_next();

    string positions_line;
    getline(input_positions_file,positions_line);

    if(!input_intensity_file.eof()) {
      if((is.x > (left_hump_end*bin_size)) &&
         (is.x < (right_hump_end*bin_size))) {
        output_good_intensity_file.write(is);
        output_good_noise_file.write(in);
        output_good_positions_file << positions_line << endl;
      } else {
        if(discarded==0) {
          output_bad_intensity_file.open(output_bad_intensityfilename);
          output_bad_noise_file    .open(output_bad_noisefilename);
          output_bad_positions_file.open(output_bad_positionsfilename);
        }
        output_bad_intensity_file.write(is);
        output_bad_noise_file.write(in);
        output_bad_positions_file << positions_line << endl;
        discarded++;
      }
    } 
  }

  out << "Discarded reads: " << discarded << endl;

  input_intensity_file.close();
  input_noise_file.close();
  
  output_good_intensity_file.close();
  output_good_noise_file.close();
  output_good_positions_file.close();

  output_bad_intensity_file.close();
  output_bad_noise_file.close();
  output_bad_positions_file.close();
}

#endif
