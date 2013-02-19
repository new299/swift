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

#ifndef SWIFT_CLUSTERFILTER_MAKEPOSITIVE_H
#define SWIFT_CLUSTERFILTER_MAKEPOSITIVE_H

#include "clusterfunctions.h"
#include <vector>
#include <iostream>
#include "ReadIntensity.h"
#include "Cluster.h"

using namespace std;

/// This class scales everything as positive, it's useful because things like purity,
/// or even the phasing correction are not calculated correctly for negative values.
template<class _prec=double>
class ClusterFilter_MakePositive {
public:

  ClusterFilter_MakePositive(const vector<Cluster<_prec> > &clusters,                        ///< Process these clusters
                             string                         source_signalid_in="RAW",        ///< Read from this signal ID
                             string                         target_signalid_in="RAW",        ///< Write to this signal ID
                             ostream                       &err_in            =std::cerr     ///< Write errors/debug info here
                            ) : source_signalid(source_signalid_in),
                                target_signalid(target_signalid_in),
                                err(err_in),
                                min_base_intensity(0),
                                add_to_base_intensity(0) {
    initialise(clusters);
  }

  /// Process a vector of clusters
  bool process(vector<Cluster<_prec> > &c) {
    for(typename vector<Cluster<_prec> >::iterator i=c.begin();i != c.end();i++) {
      process(*i);
    }
  
    return true;
  }

  /// Scale this clusters intensity in source_signalid_in as postive and write it to target_signalid_in.
  /// minimal value is calculated in constructor.
  bool process(Cluster<_prec> &c) {
    // Just add on the scaling value

    c.add_signal(target_signalid);
    c.signal(target_signalid) = c.signal(source_signalid);
    c.noise(target_signalid)  = c.noise(source_signalid);

    if(add_to_base_intensity != 0)
    for(typename Cluster<_prec>::signal_vec_type::iterator j=c.signal(target_signalid).begin();j != c.signal(target_signalid).end();j++) {
      for(int n=0;n<ReadIntensity<_prec>::base_count;n++) {
        (*j).set_base(n,(*j).get_base(n)+add_to_base_intensity);
      }
    }
  
    return true;
  }

private:
  void initialise(const vector<Cluster<_prec> > &clusters) {

    // Find smallest base value
    bool first=true;
    for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
      for(typename Cluster<_prec>::signal_vec_type::const_iterator j=(*i).const_signal(source_signalid).begin();j != (*i).const_signal(source_signalid).end();j++) {
        for(int n=0;n<ReadIntensity<_prec>::base_count;n++) {
          _prec baseval = (*j).get_base(n);
          if(baseval < min_base_intensity || first) { min_base_intensity=baseval; first=false; }
        }
      }
    }
    
    err << "Smallest base intensity value is: " << min_base_intensity << endl;
    
    if(min_base_intensity < 0) add_to_base_intensity = 0-min_base_intensity;
    else add_to_base_intensity = 0;
  }
  
  string source_signalid;                          ///< Signal ID to read data from.
  string target_signalid;                          ///< Signal ID to write data to.
  ostream &err;                                    ///< Stream to write errors/debugging info to
  
  _prec min_base_intensity;                        ///< Minimum intensity found
  _prec add_to_base_intensity;                     ///< Value to add in order to make everything positive
};

#endif
