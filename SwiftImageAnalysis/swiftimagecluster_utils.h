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

#ifndef SWIFTIMAGEANALYSIS_SWIFTIMAGECLUSTERUTILS
#define SWIFTIMAGEANALYSIS_SWIFTIMAGECLUSTERUTILS

#include <cstddef>
#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include "RLERun.h"
#include "RunLengthEncode.h"
#include <math.h>
#include "DSets.h"
#include "SwiftImageObject.h"
#include "SwiftImageCluster.h"

using namespace std;


template<class _prec>
void clear_all_offsets(vector<vector<SwiftImage<_prec> > > &images) {
  for(size_t base=0;base<images.size();base++) {
    for(size_t cycle=0;cycle<images[base].size();cycle++) {
      images[base][cycle].clear_offset();
    }
  }
}

template<class _prec>
bool cluster_isnotvalid(const vector<SwiftImageCluster<_prec> > &c) {
  
  for(typename vector<SwiftImageCluster<_prec> >::const_iterator i=c.begin();i != c.end();i++) {
    if(!(*i).isvalid()) return true;
  }
  
  return false;
}

template<class _prec>
bool cluster_mixed(const vector<SwiftImageCluster<_prec> > &c) {
  
  // if(c[0].ismixed()) return true;
  // else               return false;
  for(size_t n=0;n<c.size();n++) {
    if(c[n].ismixed()) return true;
  } 
  
  return false;

}

// remove repeated elements in a vector
template<class _t>
vector<_t> compact_vector(vector<_t> v) {
  vector<_t> result;
  
  sort(v.begin(),v.end());

  if(v.size() != 0) {
    for(unsigned int n=0;n<(v.size()-1);n++) {
      if(v[n] != v[n+1]) {
        result.push_back(v[n]);
      }
    }

    result.push_back(v[v.size()-1]);
  }
  
  return result;
}

// create a vector containing all values not in v, which are in the range 0 to <max in v>
// values in v should be positive...
template<class _t>
vector<_t> invert_vector(vector<_t> v) {
  vector<_t> result;

  v.push_back(-1);
  sort(v.begin(),v.end());

  for(unsigned int n=1;n<(v.size()-1);n++) {
    if((v[n] != (v[n+1]-1)) &&
       (v[n] != (v[n+1]  ))) {
      for(int i=v[n]+1;i<v[n+1];i++) {
        result.push_back(i);
      }
    }
  }

  return result;
}

template<class _prec>
bool validate_clusters(const vector<vector<SwiftImageCluster<_prec> > > &clusters,
                       typename SwiftImageCluster<_prec>::base_type      base) {

  for(int n=1;n<clusters.size();n++) {
    if(clusters[n].size() != clusters[n-1].size()) cerr << "ERROR: CLUSTERS DO NOT CONTAIN SAME NUMBER OF BASES" << endl;
    
    for(int na=0;na<clusters[n].size();na++) {
      if(clusters[n][na].features[base].size() == 0)     cerr << "ERROR: CLUSTER HAS NO FEATURES" << endl;
    }
  }

  return true;
}

template<class _prec>
bool validate_clusters(const vector<vector<SwiftImageCluster<_prec> > > &clusters) {

  for(size_t n=1;n<clusters.size();n++) {
    if(clusters[n].size() != clusters[n-1].size()) cerr << "ERROR: CLUSTERS DO NOT CONTAIN SAME NUMBER OF BASES" << endl;
    
//////    cerr << "Cluster " << n << endl;
    for(size_t na=0;na<clusters[n].size();na++) {
      for(size_t f=0;f<clusters[n][na].features.size();f++) {
        if(clusters[n][na].features[f].size() == 0)     cerr << "ERROR: CLUSTER HAS NO FEATURES" << endl;
      }
    }
  }

  return true;
}


/// This method joins two clusters, it should possibly be moved to SwiftImageCluster
template<class _prec>
vector<SwiftImageCluster<_prec> > join_clusters(const vector<SwiftImageCluster<_prec> > &cluster_x,
                                                const vector<SwiftImageCluster<_prec> > &cluster_y) {

  vector<SwiftImageCluster<_prec> > cluster;

  // use the size of the smaller cluster
  unsigned int size=0;
  if(cluster_x.size() > cluster_y.size()) size=cluster_y.size();
                                     else size=cluster_x.size();

  for(unsigned int n=0;n<size;n++) {
    SwiftImageCluster<_prec> c;

    for(int base=0;base < SwiftImageCluster<_prec>::base_count;base++) {
      c.features[base] = cluster_x[n].features[base]; // same as below, but faster?
    
    //c.features[base].insert(c.features[base].end(),cluster_x[n].features[base].begin(),cluster_x[n].features[base].end());
      c.features[base].insert(c.features[base].end(),cluster_y[n].features[base].begin(),cluster_y[n].features[base].end());
    }

    cluster.push_back(c);
  }

  return cluster;
}

template<class _prec>
vector<SwiftImageObject<_prec> > apply_offset(vector<SwiftImageObject<_prec> > features,
                                              SwiftImagePosition<> offset) {
  for(int n=0;n<features.size();n++) {
    features[n].apply_offset(offset);
  }

  return features;
}

template<class _prec>
vector<SwiftImageCluster<_prec> > copy_base(const vector<SwiftImageCluster<_prec> > &cluster_in,          ///< This cluster
                                            typename SwiftImageCluster<_prec>::base_type       base1,     ///< Source base
                                            typename SwiftImageCluster<_prec>::base_type       base2      ///< Destination base
                                           ) {

  vector<SwiftImageCluster<_prec> > cluster = cluster_in;

  // For all cycles
  for(typename vector<SwiftImageCluster<_prec> >::iterator i=cluster.begin();i != cluster.end();i++) {
    // Copy all features
    (*i).features[base2] = (*i).features[base1];

    for(unsigned int n=0;n<(*i).features[base2].size();n++) {
      (*i).features[base2][n].real = false;
    }
  }

  return cluster;
}

template<class _prec>
vector<SwiftImageCluster<> > replicate_existing(const    vector<SwiftImageCluster<_prec> >       &cluster,
                                                typename SwiftImageCluster<_prec>::base_type      base1,
                                                typename SwiftImageCluster<_prec>::base_type      base2) {

  if(base2 != SwiftImageCluster<_prec>::base_invalid) {
    return copy_base(cluster,base1,base2);

  } else {
    vector<SwiftImageCluster<_prec> > c = cluster;

    // Identify a base/channel with features
    int template_base=-1;
    for(int base=0;base < SwiftImageCluster<_prec>::base_count;base++) {
      if(c[0].features[base].size() != 0) {
        template_base = base;
      }
    }
    if(template_base != -1) {
      for(int base=0;base < SwiftImageCluster<_prec>::base_count;base++) {
        // Replicate existing data across all other channels
        if(c[0].features[base].size() == 0) {
          c = copy_base(c,template_base,base);
        }
      }
    } else {
      cerr << "ERROR: Cluster has no bases to replicate (not possible)" << endl;
    }
    return c;
  }
}


#endif
