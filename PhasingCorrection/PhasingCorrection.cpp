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

#include "Cluster.h"
#include "ReadIntensity.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include "clusterfilter_negativezero.h"
#include "clusterfilter_flowcell.h"
#include "clusterfilter_erode_crosstalk.h"
#include "clusterfilter_purity.h"
#include "clusterfilter_inversepurity.h"
#include "clusterfilter_offedge.h"
#include "clusterfilter_blobs.h"
#include <cmath>
#include "Timetagger.h"


//! This method generates Phasing estimates for a given cycle 
//
///! Phasing estimates are based on the amount of the called (maximum intensity) base that is incorporated
///! in to non-called (non-maximal) bases in the next (forward phasing) and previous (reverse phasing) cycles.
///! This is represented as a fraction of the base call in the current cycle.
///! Purity filtering takes place to ensure we only use good, unmixed clusters. The median phasing value across
///! pure clusters is used so that the calculation is more robust to outliers than the mean.
template<class _prec>
bool PhasingCorrection<_prec>::get_phasing(const vector<Cluster<_prec> > &clusters,int cycle,string signalid) {

  err << m_tt.str() << "Getting Phasing Estimates" << endl;
  
  // Setup filters, currently only filtering on Purity and clusters that fell off the edge of the tile,
  // but it may be worth experimenting with the other filters at some point.
  ClusterFilter_OffEdge<_prec>            filter_offedge(signalid);
  // ClusterFilter_NegativeZero<_prec>    filter_negativezero;
  // ClusterFilter_Flowcell<_prec>        filter_flowcell(clusters);
  ClusterFilter_PurityHighest<_prec>      filter_purity(clusters,12,400,signalid);
  // ClusterFilter_InversePurity<_prec>   filter_invpurity(12,0.9,signalid);
  // ClusterFilter_Blobs<_prec>           filter_blobs(clusters,5,95);
  // ClusterFilter_Erode_Crosstalk<_prec> filter_erode_ac(clusters,100,10,ReadIntensity<_prec>::base_a,ReadIntensity<_prec>::base_c);
  // ClusterFilter_Erode_Crosstalk<_prec> filter_erode_gt(clusters,100,10,ReadIntensity<_prec>::base_g,ReadIntensity<_prec>::base_t);

  int discarded=0;
  int used_bases=0;
  
  vector<vector<_prec> > all_forward_phasing(ReadIntensity<_prec>::base_count,vector<_prec>());
  vector<vector<_prec> > all_reverse_phasing(ReadIntensity<_prec>::base_count,vector<_prec>());
  for(typename ReadIntensity<_prec>::base_type b=0;b < ReadIntensity<_prec>::base_count;b++) {
    forward_phasing           [cycle][b]=0;
    forward_phasing_base_count[cycle][b]=0;
    reverse_phasing           [cycle][b]=0;
    reverse_phasing_base_count[cycle][b]=0;
  }

  // Process all clusters, get forward/reverse phasing per cluster/base add all phasing values to a vector
  for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
    
    // TODO: The logic below has been abandoned in the quest for speed. Only the purity and offedge
    // filters were in use at the time of the change. Those two filters (at least) now have an 'is_valid'
    // method which just returns a boolean, rather than making a copy of the cluster and setting its
    // valid flag, as is done by the 'process' method.
    
    // Apply filters (these mostly change the valid flag of the cluster)
    // Cluster<_prec> filtered_cluster;
    // filtered_cluster = filter_purity      .process(*i);
    // filtered_cluster = filter_offedge     .process(filtered_cluster);
    // filtered_cluster = filter_negativezero.process(filtered_cluster);
    // filtered_cluster = filter_flowcell    .process(filtered_cluster);
    // filtered_cluster = filter_invpurity   .process(filtered_cluster);
    // filtered_cluster = filter_blobs       .process(filtered_cluster);
    // filtered_cluster = filter_erode_gt    .process(filtered_cluster);
    // filtered_cluster = filter_erode_ac    .process(filtered_cluster);

    if (filter_offedge.is_valid(*i) && filter_purity.is_valid(*i)) {
      
      ReadIntensity<_prec> last_base;
      ReadIntensity<_prec> this_base;
      ReadIntensity<_prec> next_base;

      // Find next and last bases, if they exist. 
      
      bool last_base_bad = false;
      bool next_base_bad = false;

      if(cycle-1 >= 0) {
        last_base = i->const_signal(signalid)[cycle-1]; 
      } else {
        last_base_bad=true;
      }
      this_base = i->const_signal(signalid)[cycle];
      if(static_cast<unsigned int>(cycle+1) < i->const_signal(signalid).size()) {
        next_base = i->const_signal(signalid)[cycle+1];
      } else {
        next_base_bad=true;
      }

      // Forward phasing calculations
      if(!next_base_bad) {
        if(this_base.max_base() != next_base.max_base()) {
          _prec new_forwardphasing = 0;

          if(this_base.max_intensity() > 0.0001) {
            // err << "next_base.max_intensity(): " << next_base.get_base(this_base.max_base()) << " ";
            // err << "this_base.max_intensity(): " << this_base.max_intensity();
            new_forwardphasing = next_base.get_base(this_base.max_base())/this_base.max_intensity();


            all_forward_phasing[this_base.max_base()].push_back(new_forwardphasing);
            // forward_phasing           [cycle][this_base.max_base()] += new_forwardphasing;
            forward_phasing_base_count[cycle][this_base.max_base()]++;
          }
        }
      }

      // Reverse phasing calculations
      if(!last_base_bad) {
        if(this_base.max_base() != last_base.max_base()) {
          _prec new_reversephasing = 0;

          if(this_base.max_intensity() > 0.0001) {
            new_reversephasing = last_base.get_base(this_base.max_base())/this_base.max_intensity();
            all_reverse_phasing[this_base.max_base()].push_back(new_reversephasing);
            reverse_phasing_base_count[cycle][this_base.max_base()]++;
          }
        }
      }

      used_bases++;
    } else {
      discarded++;
    }

  }

  // Generate median
  // TODO: Sort this mess out (top 50 purest bases per channel?)
  err << m_tt.str() << "Sorting to generate median" << endl;
  for(typename ReadIntensity<_prec>::base_type base=0;base < ReadIntensity<_prec>::base_count;base++) {
    if (all_forward_phasing[base].size() > 0) {
      sort(all_forward_phasing[base].begin(),all_forward_phasing[base].end());
      forward_phasing[cycle][base] = all_forward_phasing[base][all_forward_phasing[base].size()/2];
    }
    if (all_reverse_phasing[base].size() > 0) {
      sort(all_reverse_phasing[base].begin(),all_reverse_phasing[base].end());
      reverse_phasing[cycle][base] = all_reverse_phasing[base][all_reverse_phasing[base].size()/2];
    }
  }


  err << m_tt.str() << "Processed Cycle: " << right << setw(2) << cycle+1 << endl;
  err << m_tt.str() << "Used Bases: " << right << setw(5) << used_bases << endl;
  err << m_tt.str() << "Discarded : " << right << setw(5) << discarded << endl;
  
  for(int n=0;n<ReadIntensity<_prec>::base_count;n++) {
    err << m_tt.str() << "Base count, first base:  " << ReadIntensity<_prec>::base_name[n] << " :"
      << right << setw(5) << forward_phasing_base_count[cycle][n] << endl;
  }
  
  for(int n=0;n<ReadIntensity<_prec>::base_count;n++) {
    err << m_tt.str() << "Base count, second base: " << ReadIntensity<_prec>::base_name[n] << " :"
      << right << setw(5) << reverse_phasing_base_count[cycle][n] << endl;
  }
  
  for(int n=0;n<ReadIntensity<_prec>::base_count;n++) {
    
    err << m_tt.str() << "Forward Phasing " << ReadIntensity<_prec>::base_name[n] << " :" << forward_phasing[cycle][n] << endl;///forward_phasing_base_count[cycle][n] << endl;
   // forward_phasing[cycle][n] = forward_phasing[cycle][n]/forward_phasing_base_count[cycle][n];
  }
  
  for(int n=0;n<ReadIntensity<_prec>::base_count;n++) {
    err << m_tt.str() << "Reverse Phasing " << ReadIntensity<_prec>::base_name[n] << " :" << reverse_phasing[cycle][n] << endl;///reverse_phasing_base_count[cycle][n] << endl;
   //  reverse_phasing[cycle][n] = reverse_phasing[cycle][n]/reverse_phasing_base_count[cycle][n];
  }

  all_phasing_zero = true;
  for(typename ReadIntensity<_prec>::base_type base=0;base < ReadIntensity<_prec>::base_count;base++) {
    if(forward_phasing[cycle][base] != 0) all_phasing_zero = false;
    if(reverse_phasing[cycle][base] != 0) all_phasing_zero = false;
  }
  
  err << m_tt.str() << "Phasing Estimation complete" << endl;
  
  return true;
}

/// Applies the current matrix to these clusters
template<class _prec>
bool PhasingCorrection<_prec>::apply_correction(vector<Cluster<_prec> > &clusters, ///< Clusters to process
                                                int cycle,                          ///< To which cycle are we applying the correction
                                                bool &thresholdmet,                 ///< Return true if we hit the threshold (and stopped)
                                                string local_source_signal,        ///< Source signal id
                                                string local_target_signal         ///< Target signal id, can be the same as source
                                               ) {
  
  err << m_tt.str() << "Applying correction" << endl;
  for(typename vector<Cluster<_prec> >::iterator i = clusters.begin();i != clusters.end();i++) {
    
    typename Cluster<_prec>::signal_vec_type old_signal;
    old_signal = (*i).const_signal(local_source_signal);
    
    typename Cluster<_prec>::signal_vec_type new_signal;

    thresholdmet=false;
    if(cycle < old_signal.size()) {
      new_signal = apply_correction(old_signal,cycle,thresholdmet);//thresholdmet is only dependant of the current phasing values for get_phasing, it should be set elsewhere.
    } else {
      new_signal = old_signal;
    }

    (*i).add_signal(local_target_signal);
    (*i).signal(local_target_signal).clear();
    (*i).signal(local_target_signal) = new_signal;

    //TODO: Inefficient
    (*i).noise(local_target_signal) = (*i).noise(local_source_signal);
  
  }
  
  err << m_tt.str() << "Correction applied" << endl;

  if(thresholdmet) return false;
  else return true;
}

/// Applies a correction to a single intensity this all need refactoring
template<class _prec>
typename Cluster<_prec>::signal_vec_type PhasingCorrection<_prec>::apply_correction(const typename Cluster<_prec>::signal_vec_type &old_signal,
                                                                                    int cycle,
                                                                                    bool &thresholdmet
                                                                                   ) {
  typename Cluster<_prec>::signal_vec_type new_signal = old_signal;
  
  // 1. Add forward phasing
  for(int b=0;b<ReadIntensity<_prec>::base_count;b++) {
    _prec use_phasing=0;
    if(forward_phasing[cycle][b] < phasing_threshold) use_phasing = forward_phasing[cycle][b];
    else  {use_phasing=phasing_threshold; thresholdmet=true;}

    if(use_phasing != 0)
    for(int na=cycle+1;(na<static_cast<int>(old_signal.size())) && (na<cycle+phasing_window);na++) {
      _prec addthis = old_signal[cycle].get_base(b)*pow(use_phasing,static_cast<int>(na-cycle));
      new_signal[cycle].set_base(b,new_signal[cycle].get_base(b) + addthis);
    }
  }
  
  // 2. Add reverse phasing
  for(int b=0;b<ReadIntensity<_prec>::base_count;b++) {
    _prec use_phasing=0;
    if(reverse_phasing[cycle][b] < phasing_threshold) use_phasing = reverse_phasing[cycle][b];
    else {use_phasing = phasing_threshold; thresholdmet=true;}

    if(use_phasing != 0)
    for(int na=cycle-1;(na>=0) && (na>=cycle-phasing_window);na--) {
      _prec addthis = old_signal[cycle].get_base(b)*pow(use_phasing,static_cast<int>(cycle-na));
      new_signal[cycle].set_base(b,new_signal[cycle].get_base(b) + addthis);
    }
  }
  
  // 3. Subtract forward phasing
  for(int b=0;b<ReadIntensity<_prec>::base_count;b++) {
    _prec use_phasing=0;
    if(forward_phasing[cycle][b] < phasing_threshold) use_phasing = forward_phasing[cycle][b];
    else {use_phasing=phasing_threshold; thresholdmet=true;}

    if(use_phasing != 0)
    for(int na=cycle+1;(na<static_cast<int>(old_signal.size())) && (na<cycle+phasing_window);na++) {
      _prec subtractthis = old_signal[cycle].get_base(b)*pow(use_phasing,static_cast<int>(na-cycle));
      new_signal[na].set_base(b,new_signal[na].get_base(b) - subtractthis);
    }
  }
  
  // 4. Subtract reverse phasing
  for(int b=0;b<ReadIntensity<_prec>::base_count;b++) {
    _prec use_phasing=0;
    if(reverse_phasing[cycle][b] < phasing_threshold) use_phasing = reverse_phasing[cycle][b];
    else {use_phasing = phasing_threshold; thresholdmet=true;}

    if(use_phasing != 0)
    for(int na=cycle-1;na>=0 && (na>=cycle-phasing_window);na--) {
      _prec subtractthis = old_signal[cycle].get_base(b)*pow(use_phasing,static_cast<int>(cycle-na));
      new_signal[na].set_base(b,new_signal[na].get_base(b) - subtractthis);
    }
  }
  
  return new_signal;
}
