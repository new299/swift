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

#ifndef SWIFT_CLUSTERFILTER_FIRSTOFFEDGE_H
#define SWIFT_CLUSTERFILTER_FIRSTOFFEDGE_H

#include <vector>
#include <iostream>
#include "Cluster.h"

using namespace std;

/// Invalidate anything that's fallen off the edge of the tile, we're not going to use it.
template<class _prec=double>
class ClusterFilter_FirstOffEdge {
public:

  ClusterFilter_FirstOffEdge(string signalid_in,int  howmany_in=1) : signalid(signalid_in), howmany(howmany_in) {
  }

  void process(vector<Cluster<_prec> > &clusters) {
    for(size_t n=0;n<clusters.size();n++) {
      process(clusters[n]);
    }
  }

  Cluster<_prec> process_create(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;

    process(co);

    return co;
  }

  void process(Cluster<_prec> &c) {

    if(c.first_off_edge(signalid,howmany)) c.valid = false;

    for(int pos=0;pos<howmany;pos++) {
      if( ((c.const_signal(signalid)[pos].get_base(ReadIntensity<_prec>::base_a) > -1) && (c.const_signal(signalid)[pos].get_base(ReadIntensity<_prec>::base_a) < 1)) ||
          ((c.const_signal(signalid)[pos].get_base(ReadIntensity<_prec>::base_c) > -1) && (c.const_signal(signalid)[pos].get_base(ReadIntensity<_prec>::base_c) < 1)) ||
          ((c.const_signal(signalid)[pos].get_base(ReadIntensity<_prec>::base_g) > -1) && (c.const_signal(signalid)[pos].get_base(ReadIntensity<_prec>::base_g) < 1)) ||
          ((c.const_signal(signalid)[pos].get_base(ReadIntensity<_prec>::base_t) > -1) && (c.const_signal(signalid)[pos].get_base(ReadIntensity<_prec>::base_t) < 1)) ) {
        c.valid=false;
      }
    }
  }

  string signalid;
  int howmany;
};

#endif
