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

#ifndef SWIFT_CLUSTERFILTER_INVERSEPURITY_H
#define SWIFT_CLUSTERFILTER_INVERSEPURITY_H

#include <string>

/// Filter based on purity, filters out reads which are very pure, marks pure clusters
/// as invalid, inpure as valid.
template<class _prec=double>
class ClusterFilter_InversePurity {
public:

  ClusterFilter_InversePurity(int num_bases_in             =12,           ///< Minimal purity with first num_bases_in bases
                              double purity_threshold_in   =0.9,          ///< Purity is less than this
                              std::string source_signal_id_in="RAW"       ///< Use this signal to determine purity
                             ) : num_bases(num_bases_in),
                                 purity_threshold(purity_threshold_in),
                                 source_signal_id(source_signal_id_in) {
  }

  /// Process a cluster, create a new cluster that is valid if inpure, otherwise invalid.
  Cluster<_prec> process(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;

    if(co.signal(source_signal_id).size() > (num_bases+1)) {
      if(co.min_purity(0,num_bases-1,source_signal_id) < purity_threshold) {
        co.valid = true;
      } else co.valid = false;
    }
    return co;
  }

  unsigned int num_bases;     ///< Number of bases over which to find minimal purity
  double purity_threshold;    ///< Threshold (purity greater than this is bad)
  string source_signal_id;    ///< Signal on which to determine purity
};

#endif
