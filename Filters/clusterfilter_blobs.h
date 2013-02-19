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

#ifndef SWIFT_CLUSTERFILTER_BLOBS_H
#define SWIFT_CLUSTERFILTER_BLOBS_H

#include <vector>
#include "ReadIntensity.h"
#include <iostream>
#include "clusterfunctions.h"

using namespace std;

/// This class is designed to filter blobs, it is quite harsh in it's filtering
/// and I wouldn't use it to filter all data, but in order to filter for use by
/// certain calculations it is useful.
/// It currently works by filtering intensity outliers, that is very weak or
/// very strong clusters, such clusters are often some from of "non-DNA" contamination.
template<class _prec=double>
class ClusterFilter_Blobs {
public:

  ClusterFilter_Blobs(const vector<Cluster<_prec> > &clusters,   ///< List of clusters in which to find upper and lower percentiles
                      int percentile_lower_limit_in = 10,        ///< Lower percentile limit
                      int percentile_upper_limit_in = 90)        ///< Upper percentile limit
                     : percentile_lower_limit(percentile_lower_limit_in),
                       percentile_upper_limit(percentile_upper_limit_in) {
    initialise(clusters);
  }

  /// Process clusters, mark those with an average maximum (call) intensity more then the upper percentile or
  /// less than the lower percentile as invalid.
  Cluster<_prec> process(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;

    _prec average = c.average_peaksignal();
    if((average > percentile_upper) || (average < percentile_lower)) {
      co.valid = false;
    } 
    
    return co;
  }
  
private:

  /// Initialisation, determines upper and lower percentile for a given set of clusters (from intensity)
  /// Q: Which signal is this using it isn't specified.
  bool initialise(const vector<Cluster<_prec> > &clusters) {
    cluster_get_percentile_peaksignal(clusters,
                                      percentile_lower_limit,
                                      percentile_upper_limit,
                                      percentile_lower,
                                      percentile_upper);
    return true;
   }

  _prec percentile_lower;           ///< The calculated value of the lower percentile limit
  _prec percentile_upper;           ///< The calculated value of the upper percentile limit
  int   percentile_lower_limit;     ///< The upper percentile (i.e. 10 for 10%)
  int   percentile_upper_limit;     ///< The lower percentile
};

#endif
