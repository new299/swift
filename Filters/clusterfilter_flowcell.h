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

#ifndef SWIFT_CLUSTERFILTER_FLOWCELL_H
#define SWIFT_CLUSTERFILTER_FLOWCELL_H

#include "clusterfunctions.h"
#include <vector>
#include <iostream>

using namespace std;

/// This class is designed to filter flowcell walls.
/// It detects a hump low to high to low transition in a series of bins
/// the transition threshold is the average*threshold, we also keep inpure
/// clusters in this region, this is a hang over from per-normalisation analysis
/// and probably screws things up now we are doing normalisation...
template<class _prec=double>
class ClusterFilter_Flowcell {
public:

  ClusterFilter_Flowcell(const vector<Cluster<_prec> > &clusters,              ///< Clusters to process
                         int                            bin_count_in=150,      ///< Number of bins to put across X-axis
                         _prec                          threshold_in=0.4,      ///< Bins count drops below/above this fraction of the average 
                         ostream                       &err_in      =std::cerr ///< Write errors here
                        ) : bin_count(bin_count_in),
                            threshold(threshold_in),
                            err(err_in) {
    initialise(clusters);
  }

  /// This processing function creates a new cluster, which is invalid if not in flowcell otherwise, unchanged.
  Cluster<_prec> process(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;
    // Now we had identified the hump end position delete all reads before the hump

    // Filter clusters based on the humps we identified.
    if((co.get_position().x > (left_hump_end*bin_size )) &&
       (co.get_position().x < (right_hump_end*bin_size))) {
      if(c.min_purity(0,11) < 0.75) {                   // Not very pure keep it
        co.valid = c.valid;
      } else {
        co.valid = false;
      }
    }
  
    return co;
  }

private:
 
  /// Create bins, based on the given set of clusters 
  void initialise(const vector<Cluster<_prec> > &clusters) {
    
    // Calculate bins
    vector<int> xbins;
    vector<int> ybins;
    bin_size=0;
    cluster_xybins(clusters,bin_count,bin_size,xbins,ybins);
  }
  
  int left_hump_end;                    ///< Left hump bin
  int right_hump_end;                   ///< Right hump bin
  int bin_count;                        ///< The number of bins to use in detection
  _prec bin_size;                       ///< Size of bin in pixels
  _prec threshold;                      ///< The threshold scaler to use in detection
  ostream &err;                         ///< Stream to write errors/debugging info to
};

#endif
