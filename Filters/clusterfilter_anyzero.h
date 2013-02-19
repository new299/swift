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

#ifndef SWIFT_CLUSTERFILTER_ANYZERO_H
#define SWIFT_CLUSTERFILTER_ANYZERO_H

#include <vector>
#include <iostream>

using namespace std;

// Invalidate anything with a 0 in the specified cycle
template<class _prec=double>
class ClusterFilter_AnyZero {
public:

  ClusterFilter_AnyZero(string signalid_in,int cycle_in=1) : signalid(signalid_in), cycle(cycle_in) {
  }

  Cluster<_prec> process(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;
    if(is_valid(c) == false) co.valid = false;
    return co;
  }

  inline bool is_valid (const Cluster<_prec> &c) {
    if(
        ((c.const_signal(signalid)[cycle].get_base(ReadIntensity<_prec>::base_a) > -1) && (c.const_signal(signalid)[0].get_base(ReadIntensity<_prec>::base_a) < 1 )) ||
        ((c.const_signal(signalid)[cycle].get_base(ReadIntensity<_prec>::base_c) > -1) && (c.const_signal(signalid)[0].get_base(ReadIntensity<_prec>::base_c) < 1 )) ||
        ((c.const_signal(signalid)[cycle].get_base(ReadIntensity<_prec>::base_g) > -1) && (c.const_signal(signalid)[0].get_base(ReadIntensity<_prec>::base_g) < 1 )) ||
        ((c.const_signal(signalid)[cycle].get_base(ReadIntensity<_prec>::base_t) > -1) && (c.const_signal(signalid)[0].get_base(ReadIntensity<_prec>::base_t) < 1 ))) {
      return false;
    } else return true;
  }

  string signalid;
  int cycle;
};

#endif
