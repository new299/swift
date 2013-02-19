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

#ifndef SWIFT_CLUSTERFILTER_ERODE_H
#define SWIFT_CLUSTERFILTER_ERODE_H

#include "clusterfunctions.h"
#include <vector>
#include <iostream>
#include "ReadIntensity.h"

using namespace std;

/// Erodes the crosstalk between bases defaults to A and C. If intensities are less than 0 this will break
/// This is a crude outlier detector for use in generating crosstalk plots. It bins A/C crosstalk and discards
/// A/C intensity pairs in regions which contain few pairs (default less than 10).
template<class _prec=double>
class ClusterFilter_Erode_Crosstalk {
public:

  ClusterFilter_Erode_Crosstalk(const vector<Cluster<_prec> >          &clusters,                                          ///< Clusters to process
                               int                                      bin_count_in = 100,                                ///< Number of bins to create across pair plot
                               int                                      threshold_in = 10,                                 ///< If this number of pairs in a bin, keep them 
                               typename ReadIntensity<_prec>::base_type base_x_in    = ReadIntensity<_prec>::base_a,       ///< First base for pair, intensity taken from first cycle of this
                               typename ReadIntensity<_prec>::base_type base_y_in    = ReadIntensity<_prec>::base_c,       ///< Second base for pair
                               string source_signalid_in                             = "RAW",                              ///< Use this signal type
                               ostream                                 &err_in       = std::cerr                           ///< Write errors here
                             ) : bin_count(bin_count_in),
                                 threshold(threshold_in),
                                 base_x(base_x_in),
                                 base_y(base_y_in),
                                 source_signalid(source_signalid_in),
                                 err(err_in) {
    initialise(clusters);
  }

  
  /// Mark oulier pairs as invalid for this cluster.
  Cluster<_prec> process(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;

    // 1. Find bin
    int bin_x = static_cast<int>((co.const_signal(source_signalid)[0].get_base(base_x)+bin_offset)/bin_size);
    int bin_y = static_cast<int>((co.const_signal(source_signalid)[0].get_base(base_y)+bin_offset)/bin_size);
    
    int friends=0;
    for(int x=-1;x<=1;x++) {
      for(int y=-1;y<=1;y++) {
        if(((bin_x+x)>0) && ((bin_x+x)<static_cast<int>(bins.size())) && ((bin_y+y)>0) && ((bin_y+y)<static_cast<int>(bins[0].size())) ) {
          if(bins[bin_x+x][bin_y+y] >= threshold) friends++;
        }
      }
    }
    if(friends <= 1) co.set_valid(false);

    return co;
  }

private:
  
  /// Create the crosstalk bins
  void initialise(const vector<Cluster<_prec> > &clusters) {
    
    // Calculate bins
    bin_size=0;
    cluster_crosstalk_bins(clusters,source_signalid,bins,base_x,base_y,bin_count,bin_offset,bin_size);
  }

  vector<vector<int> > bins;                       ///< 2D Bins containing count of reads at these crosstalk positions
  int bin_count;                                   ///< The number of bins to use in detection
  _prec bin_offset;                                ///< The number of bins to use in detection
  _prec bin_size;                                  ///< Size of bin in pixels
  int threshold;                                   ///< Minimum number of clusters per bin
  typename ReadIntensity<_prec>::base_type base_x; ///< Base in X direction to erode
  typename ReadIntensity<_prec>::base_type base_y; ///< Base in Y direction to erode
  string source_signalid;                          ///< Signal ID to read data from.
  ostream &err;                                    ///< Stream to write errors/debugging info to

};

#endif
