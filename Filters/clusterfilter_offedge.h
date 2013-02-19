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

#ifndef SWIFT_CLUSTERFILTER_OFFEDGE_H
#define SWIFT_CLUSTERFILTER_OFFEDGE_H

#include <vector>
#include <iostream>

using namespace std;

/// Invalidate anything that's fallen off the edge of the tile, we're not going to use it.
template<class _prec=double>
class ClusterFilter_OffEdge {
public:

  ClusterFilter_OffEdge(string signalid_in,int threshold_in=1) : signalid(signalid_in), threshold(threshold_in) {
  }

  void process(vector<Cluster<_prec> > &clusters) {
  
    for(size_t n=0;n<clusters.size();n++) {
      if(!is_valid(clusters[n])) clusters[n].set_valid(false);
    }
  }

  Cluster<_prec> process(const Cluster<_prec> c) {
    Cluster<_prec> co = c;
    if(is_valid(c) == false) co.valid = false;
    return co;
  }

  inline bool is_valid (const Cluster<_prec> c) {
    return ( ! c.any_off_edge(signalid,threshold));
  }

  string signalid;
  int threshold;
};

#endif
