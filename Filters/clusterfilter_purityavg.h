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

#ifndef SWIFT_CLUSTERFILTER_PURITYAVERAGE_H
#define SWIFT_CLUSTERFILTER_PURITYAVERAGE_H

#include <string>
#include <vector>

/// Filter based on purity, creates new clusters marked as valid if pure, invalid if not pure.
template<class _prec=double>
class ClusterFilter_PurityAverage {
public:

  ClusterFilter_PurityAverage(double      purity_threshold_in=0.6,   ///< Invalid if falls below this
                              size_t      purity_length_in=35,
                              std::string source_signal_id_in="RAW"  ///< Determine purity based on this signal
                             ) : purity_threshold(purity_threshold_in),
                                 purity_length(purity_length_in),
                                 source_signal_id(source_signal_id_in) {
  }

  inline void process(vector<Cluster<_prec> > &clusters) {
    for(typename vector<Cluster<_prec> >::iterator i=clusters.begin();i != clusters.end();i++) {
      process((*i));
    }
  }
  

  inline void process(Cluster<_prec> &cluster) {
    if(cluster.average_purity(source_signal_id,purity_length) > purity_threshold) {
      cluster.set_valid(true);
    } else cluster.set_valid(false);
  }

  
  /// Process cluster, creating new cluster. Valid if pure otherwise invalid.
  Cluster<_prec> process_inplace(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;

    process(c);

    return co;
  }

  double purity_threshold;   ///< Invalid if min purity falls below this threshold
  size_t purity_length;
  string source_signal_id;   ///< Determine purity based on this signal id
};

#endif
