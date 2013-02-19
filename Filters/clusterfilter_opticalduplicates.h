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

#ifndef SWIFT_CLUSTERFILTER_OPTICALDUPLICATES_H
#define SWIFT_CLUSTERFILTER_OPTICALDUPLICATES_H

#include "clusterfunctions.h"
#include <cstddef>
#include <vector>
#include <iostream>
#include "ReadIntensity.h"
#include "Cluster.h"

using namespace std;


/// This class removes what I call "Optical Duplicates" these are misidentified clusters.
/// i.e. where the image analysis as identified many clusters where there should only be
/// one, and therefore created duplicate reads.
///
/// These clusters are identified as those clusters near each other with have similar
/// sequence.
template<class _prec=double>
class ClusterFilter_OpticalDuplicates {
public:

  ClusterFilter_OpticalDuplicates(string   signalid_in            ="RAW",           ///< Detect similarity of this signal
                                  int      window_size_in         =4,               ///< Look for similar sequences in this window around each cluster
                                  int      similarity_threshold_in=2,               ///< Allow this many mismatches
                                  ostream &err_in                 =std::cerr        ///< Write errors here
                                 ) : signalid(signalid_in),
                                     window_size(window_size_in),
                                     similarity_threshold(similarity_threshold_in),
                                     err(err_in) {
  }

  /// This method processes a set of clusters. The purest cluster in a given window of a given sequence
  /// similarity this marked as valid, all else are marked as invalid. To remove optical duplicates
  /// delete all invalid clusters after processing.
  /// NOTE: This method does not mark any clusters as valid, it is assumed that all clusters are correctly
  ///       marked before processing.
  bool process(vector<Cluster<_prec> > &c) {
    
    if(c.size() == 0) return false;

    int min_x=c[0].get_position().x;
    int max_x=c[0].get_position().x;

    int min_y=c[0].get_position().y;
    int max_y=c[0].get_position().y;

    // 1. Find extreme values
    for(size_t n=0;n<c.size();n++) {
      if(c[n].get_position().x < min_x) min_x = c[n].get_position().x;
      if(c[n].get_position().x > max_x) max_x = c[n].get_position().x;
      
      if(c[n].get_position().y < min_y) min_y = c[n].get_position().y;
      if(c[n].get_position().y > max_y) max_y = c[n].get_position().y;
    }


    // 2. Build lookup image
    vector<vector<int> > lookup((max_x-min_x)+1,vector<int>((max_y-min_y)+1,-1));
    
    for(size_t n=0;n<c.size();n++) {
      lookup[c[n].get_position().x-min_x][c[n].get_position().y-min_y] = n;
    }

    // 3. Create valid/invalid vector
    vector<bool> validity(c.size(),true);

    // Iterate over clusters, if there is a purer, similar, valid
    // cluster in a X by X window, mark this cluster as invalid.

    for(size_t n=0;n<c.size();n++) {
      for(int x=c[n].get_position().x-window_size;x<=c[n].get_position().x+window_size;x++) {
        for(int y=c[n].get_position().y-window_size;y<=c[n].get_position().y+window_size;y++) {
          if((x<=max_x) && (y<=max_y) && 
             (x>=min_x) && (y>=min_y)) {
          
            int candidate = lookup[x-min_x][y-min_y];
            if((candidate != -1) && (static_cast<int>(n) != candidate)) {
              if(c[n].similarity(c[candidate],signalid) >= (static_cast<int>(c[n].signal(signalid).size())-similarity_threshold)) {
                if(c[n].min_purity(0,c[n].const_signal(signalid).size()-1,signalid) > c[candidate].min_purity(0,c[candidate].const_signal(signalid).size()-1,signalid)) validity[candidate] = false;
                                                                                       else validity[n]         = false;
              }
            }
          }
        }
      }
    }
    
    // 4. Set validity based on validity vector
    for(size_t n=0;n<c.size();n++) {
      if(validity[n] == false) c[n].valid=false;
    }

    return true;
  }

private:
  string signalid;                                 ///< Signal ID to read data from.
  int    window_size;                              ///< Size of window in which to look for similar clusters.
  int    similarity_threshold;                     ///< Clusters must be at least this similar
  ostream &err;                                    ///< Stream to write errors/debugging info to
  
};

#endif
