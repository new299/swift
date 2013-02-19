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

#ifndef SWIFT_PHASINGCORRECTION_H
#define SWIFT_PHASINGCORRECTION_H

#include "Cluster.h"
#include <iostream>
#include <vector>
#include "clusterfilter_negativezero.h"
#include "clusterfilter_purityhighest.h"
#include "Timetagger.h"
#include "clusterfunctions.h"
#include "clusterfilter_makepositive.h"

using namespace std;

template<class _prec=double>
class PhasingCorrection {
public:

  PhasingCorrection(_prec phasing_threshold_in=0.8,                  ///< Don't apply phasing greater than this (it's some kind of artifact)
                    int phasing_window_in=10,
                    string source_signalid_in="CROSSTALK_CORRECTED", ///< The signal ID to get intensities from
                    string target_signalid_in="PHASE_CORRECTED",     ///< The signal ID to write intensties to
                    ostream &err_in=std::cerr                        ///< Write debug/errors here
                   ): phasing_threshold(phasing_threshold_in),
                      source_signalid(source_signalid_in),
                      target_signalid(target_signalid_in),
                      phasing_window(phasing_window_in),
                      err(err_in) {
  }

  bool process(vector<Cluster<_prec> > &clusters) {
   
    string local_source_signalid = source_signalid;
    string local_target_signalid = target_signalid;
   
    forward_phasing.clear();
    reverse_phasing.clear();
    forward_phasing_base_count.clear();
    reverse_phasing_base_count.clear();

    forward_phasing.insert           (forward_phasing.begin()           ,(*clusters.begin()).signal(source_signalid).size(),vector<_prec>(ReadIntensity<_prec>::base_count,0));
    reverse_phasing.insert           (reverse_phasing.begin()           ,(*clusters.begin()).signal(source_signalid).size(),vector<_prec>(ReadIntensity<_prec>::base_count,0));
    forward_phasing_base_count.insert(forward_phasing_base_count.begin(),(*clusters.begin()).signal(source_signalid).size(),vector<int>  (ReadIntensity<_prec>::base_count,0));
    reverse_phasing_base_count.insert(reverse_phasing_base_count.begin(),(*clusters.begin()).signal(source_signalid).size(),vector<int>  (ReadIntensity<_prec>::base_count,0));

    // Pull out top N representative clusters
    ClusterFilter_PurityHighest<_prec> filter_purity(clusters,12,2000,source_signalid);

    vector<Cluster<_prec> > representative_clusters;
    for(typename vector<Cluster<_prec> >::iterator i=clusters.begin();i != clusters.end();i++) {
      (*i).set_valid(true);

      Cluster<_prec> c = filter_purity.process(*i);
      if(c.is_valid()) representative_clusters.push_back(c);
    }

    // We iteratively apply the correction until almost no correction was made.
    for(size_t cycle=0;cycle<(*clusters.begin()).signal(source_signalid).size();cycle++) {

      get_phasing(representative_clusters,cycle,local_source_signalid);

      // Phasing may have over corrected and made things less than 0, this will screw up subsequent cycles
      // unless we zero them.
      bool threshold_passed=false;
      if(cycle==0) {
        ClusterFilter_NegativeZero<_prec> m_clusterfilter_negativezero(local_source_signalid,local_target_signalid);
        m_clusterfilter_negativezero.process(clusters);
        m_clusterfilter_negativezero.process(representative_clusters);
      }

      if(!all_phasing_zero) {
        apply_correction(representative_clusters,cycle,threshold_passed,local_source_signalid,local_target_signalid);
        apply_correction(clusters               ,cycle,threshold_passed,local_source_signalid,local_target_signalid);
        
        local_source_signalid = target_signalid;
        local_target_signalid = target_signalid;
        // if(thresholdmet) cycle--;
        
        // Should go here or better after brace?
        // Clean up remaining zeros
        // ClusterFilter_NegativeZero<_prec> m_clusterfilter_negativezero1(local_source_signalid,local_target_signalid);
        // m_clusterfilter_negativezero1.process(clusters);
        // m_clusterfilter_negativezero1.process(representative_clusters);
      }
    }

    //ClusterFilter_NegativeZero<_prec> m_clusterfilter_negativezero(target_signalid,target_signalid);
    //m_clusterfilter_negativezero.process(clusters);

    return true;
  }

  bool get_phasing(const vector<Cluster<_prec> > &clusters,int base,string signalid);

  bool apply_correction(vector<Cluster<_prec> > &clusters,int base,bool &thresholdmet,string local_source_signalid,string local_target_signalid); ///< Applies the current matrix to these clusters
  typename Cluster<_prec>::signal_vec_type apply_correction(const typename Cluster<_prec>::signal_vec_type &r,int base,bool &thresholdmet);                            ///< Applies a correction matrix to a single intensity sequence

private:
  _prec phasing_threshold;                  ///< Don't apply phasing greater than this (artifact of some kind)

  string source_signalid;
  string target_signalid;

  vector<vector<_prec> > reverse_phasing;
  vector<vector<_prec> > forward_phasing;

  vector<vector<int>   > reverse_phasing_base_count;
  vector<vector<int>   > forward_phasing_base_count;

  int phasing_window;

  bool all_phasing_zero;

  ostream &err;                               ///< Error output will be writen here, set to cerr in constructor default
  Timetagger m_tt;

};

// This needs to be included here because this is a template class
#include "PhasingCorrection.cpp"

#endif
