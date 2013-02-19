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

#ifndef SWIFT_CLUSTERFILTER_PURITYHIGHEST_H
#define SWIFT_CLUSTERFILTER_PURITYHIGHEST_H

#include <string>
#include <algorithm>
#include "Timetagger.h"

/// Mark the top X purest clusters as valid.
template<class _prec=double>
class ClusterFilter_PurityHighest {
public:

  ClusterFilter_PurityHighest(const vector<Cluster<_prec> > &clusters,                            ///< Determine top 500 purest base on these clusters
                              int                            num_bases_in       =12,              ///< Min purity across these many initial bases
                              unsigned int                   purity_count_in    =40,              ///< How many clusters to retain after filtering
                              std::string                    source_signal_id_in="RAW",           ///< Determine purity based on this signal id.
                              ostream                       &err_in             =std::cerr        ///< Write errors here
                             ) : num_bases(num_bases_in),
                                 purity_count(purity_count_in),
                                 purity_threshold(0),
                                 source_signal_id(source_signal_id_in),
                                 err(err_in) {
 
    if(clusters.size() == 0) {
      err << "PurityHighest: There are no clusters!" << endl;
      return;
    }
    
    if(num_bases > clusters[0].const_signal(source_signal_id).size()) num_bases = clusters[0].const_signal(source_signal_id).size();

    initialise(clusters);
  }


  Cluster<_prec> process(const Cluster<_prec> &c) {
    Cluster<_prec> co = c;
    co.valid = is_valid(c);
    return co;
  }
  
  void process(vector<Cluster<_prec> > &c) {
    for(size_t n=0;n<c.size();n++) {
      c[n].set_valid(is_valid(c[n]));
    }
  }
  
  bool is_valid(const Cluster<_prec> &c) {
    
    bool ret = false;

    if(c.const_signal(source_signal_id).size() > num_bases) {
      //if(co.min_purity_greaterthaneq(0,num_bases-1,source_signal_id,purity_threshold)) {
      if(c.min_purity(0,num_bases-1,source_signal_id) >= purity_threshold) {
        ret = true;
      }
    }

    return ret;
  
  }
  
private:
  /// Determine threshold for top X purest bases  
  void initialise(const vector<Cluster<_prec> > &clusters) {
    
    //Create a vector containing all purities
    vector<_prec> purities;

    
    for(typename vector<Cluster<_prec> >::const_iterator i=clusters.begin();i != clusters.end();i++) {
      purities.push_back((*i).min_purity(0,num_bases-1,source_signal_id));
    }

    //Sort it
    sort(purities.begin(),purities.end());

    //Find entry purity_count from the end

    int position=0;
    if(purities.size() >= purity_count) position = purities.size()-purity_count;
    else {
      position = 0;
      err << m_tt.str() << "Error not enough purities: " << purities.size() << endl;
    }
    purity_threshold = purities[position];

    err << m_tt.str() << "purity_count   :  " << purity_count << endl;
    err << m_tt.str() << "purity_threshold: " << purity_threshold << endl;
  }

  unsigned int num_bases;        ///< Number of initial bases over which to determine minimum purity
  unsigned int purity_count;     ///< Number of clusters to keep?
  _prec        purity_threshold; ///< The threshold at which the above number of clusters will be kept
  string       source_signal_id; ///< The signal from which to determine purity
  ostream    &err;               ///< Write errors here
  Timetagger m_tt;

};

#endif
