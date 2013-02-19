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

#ifndef SWIFT_CLUSTERFILTER_PURITY_H
#define SWIFT_CLUSTERFILTER_PURITY_H

#include <string>

/// Filter based on purity, creates new clusters marked as valid if pure, invalid if not pure.
template<class _prec=double>
class ClusterFilter_Purity {
public:

  ClusterFilter_Purity(int         num_bases_in       =12,    ///< Min purity with this many bases
                       double      purity_threshold_in=0.6,   ///< Invalid if falls below this
                       std::string source_signal_id_in="RAW"  ///< Determine purity based on this signal
                      ) : num_bases(num_bases_in),
                          purity_threshold(purity_threshold_in),
                          source_signal_id(source_signal_id_in) {
  }

  void process(vector<Cluster<_prec> > &clusters) {
    for(size_t n=0;n<clusters.size();n++) {
      if(clusters[n].signal(source_signal_id).size() > (num_bases+1)) {
        if(clusters[n].min_purity(0,num_bases-1,source_signal_id) > purity_threshold) {
          clusters[n].valid = true;
        } else clusters[n].valid = false;
      }
    }
  }

  /// Process cluster, creating new cluster. Valid if pure otherwise invalid.
  Cluster<_prec> process(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;

    if(co.signal(source_signal_id).size() > (num_bases+1)) {
      if(co.min_purity(0,num_bases-1,source_signal_id) > purity_threshold) {
        co.valid = true;
      } else co.valid = false;
    }
    return co;
  }

  unsigned int num_bases;    ///< Min purity calculated over this many inital bases
  double purity_threshold;   ///< Invalid if min purity falls below this threshold
  string source_signal_id;   ///< Determine purity based on this signal id
};

#endif
