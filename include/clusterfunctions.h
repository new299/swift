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

#ifndef SWIFT_CLUSTERFUNCTIONS_H
#define SWIFT_CLUSTERFUNCTIONS_H

#include "Cluster.h"
#include <vector>
#include "ReadIntensity.h"
#include <algorithm>

using namespace std; // Bad me!

template<class _prec>
class Cluster_invalid
{
  public:
  
  Cluster_invalid() {
  }

  bool operator ()(const Cluster<_prec> &c) {
    bool b = (!c.is_valid());
    return b;
  }
};


template<class _prec>
void remove_invalid_clusters(vector<Cluster<_prec> > &clusters) {
  // vector<Cluster<_prec> > new_clusters;

  Cluster_invalid<_prec> pred;
  clusters.erase(remove_if(clusters.begin(),clusters.end(),pred),clusters.end());
}

/// Bins reads in the x and y direction (counting read distribution across the x and y axis)
template<class _prec>
void cluster_xybins(const vector<Cluster<_prec> > &clusters,   ///< A vector of clusters to process
                    int bin_count,                             ///< number of bins to create
                    _prec &bin_size,                           ///< returns the bin size used
                    vector<int> &xbins,                        ///< Bins in the X direction
                    vector<int> &ybins                         ///< Bins in the Y direction
                   ) {

  // find max of x and y (we'll make an equal number of bins in each direction)
  bool first=true;
  double max=0;
  for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
      if((*i).get_position().x > max || first) {max = (*i).get_position().x; first=false;}
      if((*i).get_position().y > max) max = (*i).get_position().y;
  }

  // count strips across tile
  xbins.clear();
  ybins.clear();
  xbins.insert(xbins.begin(),bin_count+1,0);
  ybins.insert(ybins.begin(),bin_count+1,0);
  bin_size = (max/(bin_count+1));
  for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
      // a little messy, but basically find the bin position (rounded) and increment
      xbins[static_cast<int>(((*i).get_position().x/(max/bin_count))+0.5)]++;
      ybins[static_cast<int>(((*i).get_position().y/(max/bin_count))+0.5)]++;
  }
}

/// Bins read in 2D space
template<class _prec>
void cluster_crosstalk_bins(const vector<Cluster<_prec> > &clusters,
                            const string &source_signalid,
                            vector<vector<int> > &bins,
                            typename ReadIntensity<double>::base_type base_x,
                            typename ReadIntensity<double>::base_type base_y,
                            int bin_count,
                            _prec &bin_offset,
                            _prec &bin_size) {

  // 1. Find max of x and y
  // find max of x and y (we'll make an equal size bins in each direction)
  bool firstmin=true;
  bool firstmax=true;
  double min=0;
  double max=0;
  for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
      if((*i).const_signal(source_signalid)[0].get_base(base_x) > max || firstmax) {max = (*i).const_signal(source_signalid)[0].get_base(base_x); firstmax=false;}
      if((*i).const_signal(source_signalid)[0].get_base(base_y) > max            )  max = (*i).const_signal(source_signalid)[0].get_base(base_y);
      
      if((*i).const_signal(source_signalid)[0].get_base(base_x) < min || firstmin) {min = (*i).const_signal(source_signalid)[0].get_base(base_x); firstmin=false;}
      if((*i).const_signal(source_signalid)[0].get_base(base_y) < min            )  min = (*i).const_signal(source_signalid)[0].get_base(base_y);
  }

  if(min<0) bin_offset= 0-min;
  else bin_offset=0;

  bins.clear();
 
  // Make 2D bins
  for(int n=0;n<bin_count;n++) {
    vector<int> v(bin_count,0);
    bins.push_back(v);
  }

  bin_size = (bin_offset+max)/(bin_count-1);
  for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
    
    _prec xbase = (*i).const_signal(source_signalid)[0].get_base(base_x);
    _prec ybase = (*i).const_signal(source_signalid)[0].get_base(base_y);

    if((!isnan(xbase)) && (!isnan(ybase))) {
      int x_bin = static_cast<int>((xbase+bin_offset)/bin_size);
      int y_bin = static_cast<int>((ybase+bin_offset)/bin_size);

      bins[x_bin][y_bin]++;
    }
  }
}

/// Locates a "hump" in read distribution across an axis (uses bins generated above)
template<class _prec>
int cluster_find_hump(const vector<int> &bins, ///< A vector containing the number of items in each bin
                      _prec threshold,         ///< Scaler for average to create threshold for high/low transitions
                      int start,               ///< Start processing from this bin
                      int end,                 ///< Stop at this bin
                      int direction,           ///< Increment in this direction i.e. 1 = 0->1, -1 = 1->0
                      ostream &out             ///< Write errors/debugging info to this stream
                     ) {

  // find max
  _prec average_in_bins=0;
  for(vector<int>::const_iterator i = bins.begin();i != bins.end();i++) {
    average_in_bins += (*i);
  }
  average_in_bins = average_in_bins/bins.size();
  out << "average: " << average_in_bins << endl;

  bool low=false;
  bool high=false;
  bool first=true;

  bool high_low_trans=false;
  bool low_high_trans=false;

  // Walk across bins, detect transition
  int hump_end=start;
  for(int n=start;(n != end) && (hump_end == 0);n=n+direction) {
    if(bins[n] < (average_in_bins*threshold)) {
      low=true;
      if(!first) {
        if(high) {high_low_trans = true;}
        high=false;
      }
      first=false;
    }
    else {
      high=true;
      if(!first) {
        if(low) {low_high_trans = true;}
        low=false;
      }
      first=false;
   }
   if(low_high_trans && high_low_trans) {
     hump_end=n;
   }
  }

  return hump_end;
}

template<class _prec>
bool cluster_get_percentile(const vector<Cluster<_prec> > &clusters,        ///< Clusters to use
                            typename ReadIntensity<_prec>::base_type base,  ///< Get the percentile limit of intensities for this base
                            int percentile_lower_limit,
                            int percentile_upper_limit,
                            double &percentile_lower,
                            double &percentile_upper
                           ) {
  // Pull out all the intensities for this base

  vector<_prec> base_intensities;

  for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
    base_intensities.push_back((*i).signal().get_base(base));
  }

  sort(base_intensities.begin(),base_intensities.end());

  percentile_lower = base_intensities[base_intensities.size()*(static_cast<double>(percentile_lower_limit)/100)];
  percentile_upper = base_intensities[base_intensities.size()*(static_cast<double>(percentile_upper_limit)/100)];

  return true;
}

template<class _prec>
bool cluster_get_percentile_peaksignal(const vector<Cluster<_prec> > &clusters,        ///< Clusters to use
                                       int percentile_lower_limit,
                                       int percentile_upper_limit,
                                       double &percentile_lower,
                                       double &percentile_upper
                                      ) {
  // Pull out all the intensities for this base

  vector<_prec> base_intensities;

  for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
    base_intensities.push_back((*i).average_peaksignal());
  }

  sort(base_intensities.begin(),base_intensities.end());

  percentile_lower = base_intensities[static_cast<int>(base_intensities.size()*(static_cast<double>(percentile_lower_limit)/100))];
  percentile_upper = base_intensities[static_cast<int>(base_intensities.size()*(static_cast<double>(percentile_upper_limit)/100))];

  return true;
}

template<class _prec>
inline void clear_cluster_signal(vector<Cluster<_prec> > &clusters,const string &signal_id) {
  for(typename vector<Cluster<_prec> >::iterator i=clusters.begin();i != clusters.end();i++) {
    (*i).delete_signal(signal_id);
  }
}

template<class _prec>
inline void clear_cluster_validity(vector<Cluster<_prec> > &clusters, bool validity) {
  for(typename vector<Cluster<_prec> >::iterator i=clusters.begin();i != clusters.end();i++) {
    (*i).set_valid(validity);
  }
}


#endif
