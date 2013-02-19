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

#ifndef INTENSITYTOOLS_INTENSITYFILTER_FLOWCELL_H
#define INTENSITYTOOLS_INTENSITYFITLER_FLOWCELL_H

#include "intensity2xybins.h"
#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;

/// This class is designed to filter flowcell walls.
/// It detects a hump low to high to low transition in a series of bins
/// the transition threshold is the average*threshold
class IntensityFilter_Flowcell {
public:

  IntensityFilter_Flowcell(IntensityReader &intensity_file_in,
                          int bin_count_in,
                          double threshold_in,
                          ostream &err_in
                         ) : input_intensity_file(&intensity_file_in),
                             bin_count(bin_count_in),
                             threshold(threshold_in),
                             err(err_in) {
    initialise();
  }

  void initialise() {
    input_intensity_file->reopen();
    
    // Calculate bins
    vector<int> xbins;
    vector<int> ybins;
    bin_size=0;
    intensity2xybins(*input_intensity_file,bin_count,bin_size,xbins,ybins);

    // find hump in xbins, a hump means the bins rise above 25% of the maximum value and then fall below 25% again within the first 25% of bins.
    
    left_hump_end  = find_hump(xbins,threshold,0,xbins.size()/4,1,err);
    right_hump_end = find_hump(xbins,threshold,xbins.size(),xbins.size()-(xbins.size()/4),-1,err);
   
    // Push these out a little to avoid numerical errors removing some reads.
    if(left_hump_end == 0) left_hump_end = -1;
    if(right_hump_end == bin_count) right_hump_end = bin_count+1;
    
    err << "Using Left Hump End Bin: " << left_hump_end << endl;
    err << "Using Left Hump End Pos: " << left_hump_end*bin_size << endl;
    
    err << "Using Right Hump End Bin: " << right_hump_end << endl;
    err << "Using Right Hump End Pos: " << right_hump_end*bin_size << endl;
  }

  IntensitySequence filter(IntensitySequence is) {
    // Now we had identified the hump end position delete all reads before the hump

    // Filter clusters based on the humps we identified.
    int discarded=0;
    if(((is.x > (left_hump_end*bin_size)) &&
       (is.x < (right_hump_end*bin_size))) ||
       (is.min_purity(0,11) < 0.75)) {                   // Not very pure keep it
      is.valid = true;
    } else {
      is.valid = false;
    }
  
    return is;
  }

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

private:
  int left_hump_end;                     ///< Left hump bin
  int right_hump_end;                    ///< Right hump bin
  int bin_count;                         ///< The number of bins to use in detection
  double bin_size;                       ///< Size of bin in pixels
  double threshold;                      ///< The threshold scaler to use in detection
  ostream &err;                          ///< Stream to write errors/debugging info to
  IntensityReader *input_intensity_file; ///< File which we are reading from
};

#endif
