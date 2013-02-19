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

#ifndef SWIFT_BASECALLER_H
#define SWIFT_BASECALLER_H

#include "Cluster.h"
#include <iostream>
#include <vector>
#include "ScoredSequence.h"

using namespace std;

template<class _prec=double>
class BaseCaller {
public:

  BaseCaller(string source_signalid_in  ="PHASE_CORRECTED", ///< The signal ID to get intensities from
             string target_sequenceid_in="BASECALL",        ///< The sequence ID to write basecalls to
             ostream &err_in=std::cerr                      ///< Write debug/errors here
            ): source_signalid(source_signalid_in),
               target_sequenceid(target_sequenceid_in),
               err(err_in) {
  }

  bool process(vector<Cluster<_prec> > &clusters) {
    for(typename vector<Cluster<_prec> >::iterator i = clusters.begin();i != clusters.end();i++) {
      
      (*i).add_sequence(target_sequenceid);// if it doesn't exist add it

      (*i).sequence(target_sequenceid).sequence().clear();
      (*i).sequence(target_sequenceid).quality() .clear();
      bool last_offedge=false;
      for(typename Cluster<_prec>::signal_vec_type::const_iterator j=(*i).const_signal(source_signalid).begin();j != (*i).const_signal(source_signalid).end();j++) {

        // More than one base off the edge calls an N otherwise call max peak.
        // Why more than one? because phasing can compensate to an extent
        // TODO: Why do I have to use this messy cast here? remove it
        if((*j).off_edge() && last_offedge) {
          (*i).sequence(target_sequenceid).sequence().push_back(static_cast<int>(Cluster<_prec>::sequence_type::base_n)); 
          (*i).sequence(target_sequenceid).quality().push_back(1); 
        } else {
          // TODO: The following line is bad, it relies of ReadIntensity and ScoredSequence seting base consts the same
          (*i).sequence(target_sequenceid).sequence().push_back((*j).max_base()); 
          (*i).sequence(target_sequenceid).quality().push_back((*j).get_quality()); 
          last_offedge=false;
        }
        if((*j).off_edge()) last_offedge=true;
        
      }
    }
    
    return true;
  }

private:

  string source_signalid;
  string target_sequenceid;

  ostream &err;                               ///< Error output will be writen here, set to cerr in constructor default
};

// This needs to be included here because this is a template class
// #include "BaseCaller.cpp"

#endif
