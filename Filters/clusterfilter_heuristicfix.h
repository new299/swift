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

#ifndef SWIFT_CLUSTERFILTER_HEURISTICFIX_H
#define SWIFT_CLUSTERFILTER_HEURISTICFIX_H

#include <vector>
#include <iostream>
#include "ReadIntensity.h"
#include "Cluster.h"
#include "clusterfilter_purity.h"

using namespace std;

/// This class is designed to perform a number of heuristic fixes to basecalls.
/// Currently it does very little. The only changed that has worked (and only
/// occassionally) is to reset the basecall if the phasing correction altered the
/// call in the first 12 cycles. This isn't used in the current Swift driver.
template<class _prec=double>
class ClusterFilter_HeuristicFix {
public:

  ClusterFilter_HeuristicFix(string    source_signalid_in="CROSSTALK_CORRECTED", ///< Read signals from this ID
                             string    target_signalid_in="NORMALISE_CORRECTED", ///< Write signals to this ID
                              ostream &err_in            =std::cerr              ///< Write debug/errors here.
                             ) : source_signalid(source_signalid_in),
                                 target_signalid(target_signalid_in),
                                 err(err_in) {
  }

  /// Process a whole vector of clusters
  bool process(vector<Cluster<_prec> > &c) {
    for(typename vector<Cluster<_prec> >::iterator i=c.begin();i != c.end();i++) {
      process(*i);
    }
  
    return true;
  }

  /// Processes a cluster, read signal from source_signalid, and write new signal to target_signalid
  bool process(Cluster<_prec> &c) {

    c.add_signal(target_signalid);
    c.signal(target_signalid) = c.signal(source_signalid);
    c.noise(target_signalid)  = c.noise(source_signalid);

    int max_cycles=12;
    if(static_cast<unsigned int>(max_cycles) > c.const_signal(source_signalid).size()) max_cycles = c.const_signal(source_signalid).size();
/*    for(int cycle=0;cycle < max_cycles;cycle++) {
      if(c.signal(target_signalid)[cycle].maxdiff() < 100) {
        if(((c.signal(target_signalid)[cycle].    max_base() == ReadIntensity<_prec>::base_a) &&
            (c.signal(target_signalid)[cycle].sub_max_base() == ReadIntensity<_prec>::base_t)) ||
           ((c.signal(target_signalid)[cycle].    max_base() == ReadIntensity<_prec>::base_g) &&
            (c.signal(target_signalid)[cycle].sub_max_base() == ReadIntensity<_prec>::base_c)))
        c.signal(target_signalid)[cycle] = c.signal("RAW")[cycle];
      }
    }
*/    
    // works a little, sometimes not at all...
    if(c.min_purity(0,11,source_signalid) < 0.56)
    for(int cycle=0;cycle < max_cycles;cycle++) {
      if(c.signal(source_signalid)[cycle].purity() < 0.56)
      if(c.signal("PHASE_CORRECTED")[cycle].    max_base() != c.signal("MAKEPOSITIVE_CORRECTED")[cycle].    max_base()) {
        c.signal(target_signalid)[cycle] = c.signal("MAKEPOSITIVE_CORRECTED")[cycle];
      }
    }
    
    // doesn't work
    /*
    for(unsigned int cycle=0;cycle < c.const_signal(source_signalid).size();cycle ++) {
      if(c.signal("MAKEPOSITIVE_CORRECTED")[cycle].purity() > 0.6) {
        c.signal(target_signalid)[cycle] = c.signal("MAKEPOSITIVE_CORRECTED")[cycle];
      }
    }*/
  
//  needs trying again  , still didn't work...
/*    if(c.min_purity(0,11) < 0.55) {
      c.signal(target_signalid) = c.signal("MAKEPOSITIVE_CORRECTED");
    }*/
    
    return true;
  }

private:
  string source_signalid;                        ///< Signal ID to read data from.
  string target_signalid;                        ///< Signal ID to write data to.
  ostream &err;                                  ///< Stream to write errors/debugging info to
  
};

#endif
