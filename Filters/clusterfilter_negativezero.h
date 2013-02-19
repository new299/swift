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

#ifndef SWIFT_CLUSTERFILTER_NEGATIVEZERO_H
#define SWIFT_CLUSTERFILTER_NEGATIVEZERO_H

#include "clusterfunctions.h"
#include <vector>
#include <iostream>
#include "ReadIntensity.h"
#include "Cluster.h"

using namespace std;

/// Any negative value is reset to zero. This is useful to discard overcorrection.
/// Negative values cause problems for things like purity, and phasing corrections.
/// So if you've overcorrected and created a negative value, set this value to zero.
template<class _prec=double>
class ClusterFilter_NegativeZero {
public:

  ClusterFilter_NegativeZero(string source_signalid_in  ="RAW",            ///< Read signal from here
                             string target_signalid_in  ="RAW",            ///< Write signal here
                             ostream &err_in            =std::cerr         ///< Write errors here
                            ) : source_signalid(source_signalid_in),
                                target_signalid(target_signalid_in),
                                err(err_in) {
  }

  /// Process a vector of clusters
  bool process(vector<Cluster<_prec> > &c) {
    for(typename vector<Cluster<_prec> >::iterator i=c.begin();i != c.end();i++) {
      process(*i);
    }
  
    return true;
  }

  /// Process a single cluster (set negative source_signalid_in values to zero and place in target_signalid_in)
  bool process(Cluster<_prec> &c) {

    c.add_signal(target_signalid);
    c.signal(target_signalid) = c.signal(source_signalid);
    c.noise(target_signalid)  = c.noise(source_signalid);

    for(typename Cluster<_prec>::signal_vec_type::iterator j=c.signal(target_signalid).begin();j != c.signal(target_signalid).end();j++) {
      for(int n=0;n<ReadIntensity<_prec>::base_count;n++) {
        if((*j).get_base(n) < 0)  (*j).set_base(n,0);
      }
    }
  
    return true;
  }

private:
  string source_signalid;                          ///< Signal ID to read data from.
  string target_signalid;                          ///< Signal ID to write data to.
  ostream &err;                                    ///< Stream to write errors/debugging info to
  
};

#endif
