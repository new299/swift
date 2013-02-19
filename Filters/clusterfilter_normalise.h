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

#ifndef SWIFT_CLUSTERFILTER_NORMALISE_H
#define SWIFT_CLUSTERFILTER_NORMALISE_H

#include <vector>
#include <iostream>
#include <iomanip>
#include "ReadIntensity.h"
#include "Cluster.h"
#include "Timetagger.h"

using namespace std;

/// Normalise signals, currently against median
template<class _prec=double>
class ClusterFilter_Normalise {
public:

  ClusterFilter_Normalise(const vector<Cluster<_prec> > &clusters,                                  ///< Determine median from here
                          string                         source_signalid_in="CROSSTALK_CORRECTED",  ///< Read signals from here
                          string                         target_signalid_in="NORMALISE_CORRECTED",  ///< Write signals here
                          ostream                       &err_in            =std::cerr               ///< Write errors here
                         ) : source_signalid(source_signalid_in),
                             target_signalid(target_signalid_in),
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

  /// Process a single cluster, normalisation finds the median for the channel and subtracts this from everything, meaning
  /// that intensities are now measured as the deviation from the median. This seems to compensate for a given background
  /// illuminate present in a channel. I've tried dividing by the median, to compensate for illumation range differences,
  /// however this didn't help. I should try dividing by the SD.
  bool process(Cluster<_prec> &c) {
    // Just add on the scaling value

    c.add_signal(target_signalid);
    c.signal(target_signalid) = c.signal(source_signalid);
    c.noise(target_signalid)  = c.noise(source_signalid);

    for(unsigned int cycle=0;cycle < c.const_signal(source_signalid).size();cycle++) {
      for(int base=0;base<ReadIntensity<_prec>::base_count;base++) {
        c.signal(target_signalid)[cycle].set_base(base,c.signal(source_signalid)[cycle].get_base(base)-average_base_intensity[cycle][base]);
      }
    }
    
    return true;
  }

private:
  /// Find the median intensity for each cycle
  void initialise(const vector<Cluster<_prec> > &clusters) {

    average_base_intensity.clear();

    // Add all clusters, all clusters MUST have the same number of cycles.
    // TODO: a vector of clusters should probably be its own object, which
    // TODO: has a method to return the number maximum number of cycles in
    // TODO: all clusters
   
    average_base_intensity.insert(average_base_intensity.begin(),clusters[0].const_signal(source_signalid).size(),vector<_prec>(ReadIntensity<_prec>::base_count,0));

    // Find average base value per cycle
    size_t cycle_count = (*(clusters.begin())).const_signal(source_signalid).size();

    for(int base=0;base<ReadIntensity<_prec>::base_count;base++) {
      for(size_t cycle=0;cycle < cycle_count;cycle++) {   // need a better way of getting cycle count, max cycles
        
        vector<_prec> all_average_base_intensity;
        for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {

          if(cycle < (*i).const_signal(source_signalid).size()) {
          
            //Ignore negative values!!!
            if(!((*i).const_signal(source_signalid)[cycle].get_base(base) <= 0)) {
              all_average_base_intensity.push_back((*i).const_signal(source_signalid)[cycle].get_base(base));
            }
          }

        }

        if(cycle < average_base_intensity.size()) {
          sort(all_average_base_intensity.begin(),all_average_base_intensity.end());
         
          if((all_average_base_intensity.size()/2) < all_average_base_intensity.size()) {
            average_base_intensity[cycle][base] = all_average_base_intensity[ all_average_base_intensity.size()/2 ];
            cout << m_tt.str() << "median cycle: " << right << setw(2) << cycle+1 
              << " base " << ReadIntensity<>::base_name[base] << ": " << average_base_intensity[cycle][base] << endl;
          }
        }
      }
    }

    // Calculate average, now using median...
    //for(typename vector<vector<_prec> >::iterator i = average_base_intensity.begin();i != average_base_intensity.end();i++) {
    //  for(typename vector<_prec>::iterator j = (*i).begin();j != (*i).end();j++) {
    //    (*j) = (*j)/clusters.size();
    //  }
    //}
  }
  
  string source_signalid;                        ///< Signal ID to read data from.
  string target_signalid;                        ///< Signal ID to write data to.
  ostream &err;                                  ///< Stream to write errors/debugging info to
  Timetagger m_tt;

  vector<vector<_prec> > average_base_intensity; ///< Average base intensity cycle/base
  
};

#endif
