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

#ifndef SWIFT_CLUSTERFILTER_SCALEQUALITY_H
#define SWIFT_CLUSTERFILTER_SCALEQUALITY_H

#include <vector>
#include <iostream>
#include "ReadIntensity.h"
#include "Cluster.h"

using namespace std;


/// Rescale quality probabilities over the range 0->1
template<class _prec=double>
class ClusterFilter_ScaleQuality {
public:

  ClusterFilter_ScaleQuality(const vector<Cluster<_prec> > &clusters, ///< Process these clusters
                             string source_sequenceid_in="BASECALLF", ///< Read scored sequences from here
                             string target_sequenceid_in="BASECALLF", ///< Write scored sequences to here
                             ostream &err_in=std::cerr                ///< Write errors here
                         ) : source_sequenceid(source_sequenceid_in),
                             target_sequenceid(target_sequenceid_in),
                             err(err_in) {
    initialise(clusters);
  }

  /// Process a vector of clusters
  bool process(vector<Cluster<_prec> > &c) {
    for(typename vector<Cluster<_prec> >::iterator i=c.begin();i != c.end();i++) {
      process(*i);
    }
  
    return true;
  }

  /// Process a single cluster (read from source_sequenceid_in write to target_sequenceid_in).
  bool process(Cluster<_prec> &c) {

    c.add_sequence(target_sequenceid);
    c.sequence(target_sequenceid) = c.sequence(source_sequenceid);

    for(unsigned int cycle=0;cycle < c.sequence(source_sequenceid).quality().size();cycle++) {
      c.sequence(target_sequenceid).quality()[cycle] = c.sequence(source_sequenceid).quality()[cycle]/max_quality;
    }
    
    return true;
  }

private:
  
  void initialise(const vector<Cluster<_prec> > &clusters) {

    max_quality = 0;
    // Add all clusters, all clusters MUST have the same number of cycles.
    // TODO: a vector of clusters should probably be it's own object, which
    // TODO: has a method to return the number maximum number of cycles in
    // TODO: all clusters
   
    // Find average base value per cycle
    bool first=true;
    for(unsigned int cycle=0;cycle < (*(clusters.begin())).const_sequence(source_sequenceid).const_quality().size();cycle++) {   // need a better way of getting cycle count, max cycles
        
      for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
        _prec new_quality = (*i).const_sequence(source_sequenceid).const_quality()[cycle];
        if((new_quality > max_quality) || first) {
          max_quality = new_quality;
          first = false;
        }
      }
    }

    err << "Max quality: " << max_quality << endl;
  }

  string source_sequenceid;                      ///< Signal ID to read data from.
  string target_sequenceid;                      ///< Signal ID to write data to.
  ostream &err;                                  ///< Stream to write errors/debugging info to
  _prec max_quality;                             ///< maximum quality value

};

#endif
