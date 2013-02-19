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

#ifndef SWIFT_PRBBASECALLER_H
#define SWIFT_PRBBASECALLER_H

#include "Cluster.h"
#include <iostream>
#include <vector>
#include "ProbabilitySequence.h"

using namespace std;

template<class _prec=float,class _prbprec=float>
class PrbBaseCaller {
public:

  PrbBaseCaller(string source_signalid_in  ="PHASE_CORRECTED", ///< The signal ID to get intensities from
                string target_sequenceid_in="BASECALL",        ///< The sequence ID to write basecalls to
                ostream &err_in=std::cerr                      ///< Write debug/errors here
               ): source_signalid(source_signalid_in),
                  target_sequenceid(target_sequenceid_in),
                  err(err_in) {
  }

  bool process(Cluster<_prec> &cluster,
               bool delete_signal) {
    cluster.add_sequence(target_sequenceid);// if it doesn't exist add it

    cluster.sequence(target_sequenceid).sequence().clear();
    bool last_offedge=false;
    for(typename Cluster<_prec>::signal_vec_type::const_iterator j=cluster.const_signal(source_signalid).begin();j != cluster.const_signal(source_signalid).end();j++) {

      // More than one base off the edge calls an N otherwise call max peak.
      // Why more than one? because phasing can compensate to an extent
      // TODO: Why do I have to use this messy cast here? remove it
      if((*j).off_edge() && last_offedge) {
        // hum....
        BaseProbability<_prbprec> prb;
        for(size_t n=0;n < BaseProbability<_prbprec>::base_count;n++) prb[n] = 1/static_cast<_prbprec>(ProbabilitySequence<_prbprec>::base_count);

        cluster.sequence(target_sequenceid).sequence().push_back(prb);
      } else {
        // TODO: The following line is bad, it relies of ReadIntensity and ProbabilitySequence seting base consts the same
        
        _prec intensity_sum=0;
        for(size_t n=0;n < ProbabilitySequence<_prec>::base_count;n++) {
          if((*j).get_base(n) > 0) intensity_sum += (*j).get_base(n);
        }

        BaseProbability<_prbprec> prb;
        for(size_t n=0;n < ProbabilitySequence<_prbprec>::base_count;n++) {
           if(intensity_sum > 0) {
             if((*j).get_base(n) > 0) prb[n] = (*j).get_base(n)/intensity_sum;
                                 else prb[n] = 0;
           } else prb[n] = 0;
        }
       
        cluster.sequence(target_sequenceid).sequence().push_back(prb);

        last_offedge=false;
      }
      if((*j).off_edge()) last_offedge=true;
    }

    if(delete_signal) {
      cluster.delete_signal(source_signalid);
    }

    return true;
  }

  // Process a set of clusters
  bool process(vector<Cluster<_prec> > &clusters,
               bool delete_signal=false 
              ) {
    for(typename vector<Cluster<_prec> >::iterator i = clusters.begin();i != clusters.end();i++) {
      process((*i),delete_signal);
    }
    
    return true;
  }

  vector<Cluster<_prec> > process_and_reallocate(vector<Cluster<_prec> > &clusters) {
    vector<Cluster<_prec> > new_clusters;
    
    for(int n = clusters.size()-1;n >= 0;n--) {
      process(clusters[n],true);

      new_clusters.push_back(clusters[n]);
      clusters.erase(clusters.begin()+n);
    }

    return new_clusters;
  }

private:

  string source_signalid;
  string target_sequenceid;

  ostream &err;                               ///< Error output will be writen here, set to cerr in constructor default
};

// This needs to be included here because this is a template class
// #include "BaseCaller.cpp"

#endif
