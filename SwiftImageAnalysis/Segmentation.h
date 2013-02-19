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

#ifndef SWIFTIMAGEANALYSIS_SEGMENTATION
#define SWIFTIMAGEANALYSIS_SEGMENTATION

#include <cstddef>
#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include "RLERun.h"
#include "RunLengthEncode.h"
#include <math.h>
#include "DSets.h"
#include "SwiftImageObject.h"
#include "Timetagger.h"
#include "EuclideanDistanceMap.h"
#include "Watershed.h"
#include "Invert.h"

using namespace std;

/// This class segments a binary image returning a bunch of SwiftImageObjects
template <class _prec=uint16,class _threshold_prec=uint16,class _image_cluster_prec=double>
class Segmentation {
public:
  
  Segmentation(bool use_watershed_in=false, ostream &err_in=std::cerr) : use_watershed(use_watershed_in), err(err_in) {
  }

  SwiftImage<int> rle_lookup(const vector<RLERun<> > &runs,int image_width,int image_height) {
    
    // TODO: These offsets (+200s) should be based on maximum offset trials in image analysis
    SwiftImage<int> lookup(image_width+200,image_height+200,-1);
    lookup.apply_offset(SwiftImagePosition<>(100,100));

    
    for(size_t n=0;n<runs.size();n++) {
      int y = runs[n].pos.y;
      for(int x=runs[n].pos.x; x<(runs[n].pos.x+runs[n].length) ;x++) {
        lookup(x,y) = n;
      }
    }

    return lookup;
  }

  vector<SwiftImageCluster<_image_cluster_prec> > process(const vector<vector<SwiftImage<_threshold_prec> > > &thresholded,
                                       const vector<vector<SwiftImage<_prec> > >         &images,
                                       int max_cycles=6) {
    
    SwiftImage<int> lookup(images[0][0].image_width()+200,images[0][0].image_height()+200,-1);
    lookup.apply_offset(SwiftImagePosition<>(100,100));
    
    vector<SwiftImageCluster<_image_cluster_prec> > clusters;

    if(max_cycles==-1) max_cycles = thresholded[0].size();

    for(int base=0;base<4;base++) {
      for(size_t cycle=0;(cycle<thresholded[base].size()) && (cycle < static_cast<size_t>(max_cycles));cycle++) {
        err << m_tt.str() << "Segmenting: " << base_list[base] << " " << cycle+1 << endl;

        // 4. Segment images
        
        // Watershed
        SwiftImage<_threshold_prec> i6 = thresholded[base][cycle];
        
        if(use_watershed) {
          EuclideanDistanceMap<_threshold_prec> edm;
          Watershed<_threshold_prec> wat;
          Invert<_threshold_prec> inv;

          SwiftImage<_threshold_prec> i3 = edm.process(thresholded[base][cycle]);
          SwiftImage<_threshold_prec> i4 = inv.process(i3);
          SwiftImage<_threshold_prec> i5 = wat.process(i4);
        } 

        vector<SwiftImageCluster<_image_cluster_prec> > cur_clusters = process(i6,images,lookup);

        clusters.insert(clusters.begin(),cur_clusters.begin(),cur_clusters.end());
        
        err << m_tt.str() << "Cycle " << cycle+1 << " Segmentation complete: " << cur_clusters.size() << " clusters" << endl;
      }
    }

    return clusters;
  }
 
  template<class _prec1>
  vector<SwiftImageCluster<_image_cluster_prec> > process(const SwiftImage<_prec1> &source,
                                       const vector<vector<SwiftImage<_prec> > > &images,
                                       SwiftImage<int> &lookup_c
                                      ) {
   
    RunLengthEncode<_prec1> rle;

    vector<RLERun<> > runs = rle.process(source);
    DSets sets(runs.size());

    // Generate lookup image, -1 means there is no run here
    SwiftImage<int> lookup = rle_lookup(runs,source.image_width(),source.image_height());

    // Iterate over runs, find adjacent runs in lookup image and join them
    for(vector<RLERun<> >::iterator i = runs.begin();i != runs.end();i++) {
      
      int y = (*i).pos.y;
      for(int x=(*i).pos.x;x < (*i).pos.x+(*i).length;x++) {
       
        if(y-1 >= 0) {
          if(lookup(x,y-1) >= 0) {
            // combine these features
            if(lookup(x,y-1) != lookup(x,y)) {
              sets.makeparent(lookup(x,y),lookup(x,y-1));
              x=runs[lookup(x,y-1)].pos.x+runs[lookup(x,y-1)].length;
            }
          }
        }
      }
    }

    // Create a bunch of features

    vector<SwiftImageObject<_image_cluster_prec> > features(runs.size(),SwiftImageObject<_image_cluster_prec>());

    for(unsigned int n=0;n<runs.size();n++) {
      int c = sets.getparent(n); // getparent returns the canonical element for this set

      features[c].pixels.push_back(runs[n]);
    }

    // Remove all empty features
    vector<SwiftImageObject<_image_cluster_prec> > features_real;
    for(typename vector<SwiftImageObject<_image_cluster_prec> >::iterator i=features.begin();i != features.end();i++) {
      if((*i).pixels.size() != 0) features_real.push_back(*i);
    }
   
    // Create clusters
    //TODO: the logic here isn't entirely correct... we could end up not adding a new cluster, but removing an existing one

    vector<SwiftImageCluster<_image_cluster_prec> > clusters;
    for(typename vector<SwiftImageObject<_image_cluster_prec> >::iterator i=features_real.begin();i != features_real.end();i++) {
      
      vector<int> overlaps = (*i).find_in_image(lookup_c);

      if(overlaps.size() == 0) {
        (*i).set_image(lookup_c,1);

        SwiftImageCluster<_image_cluster_prec> c(*i);
        clusters.push_back(c);
      }
    }
   
    return clusters; 
  }

private:

  bool use_watershed;
  ostream &err;
  Timetagger m_tt;
  static const char base_list[];

};   // end of class definition

template <class _prec,class _threshold_prec,class _image_cluster_prec>
const char Segmentation<_prec,_threshold_prec,_image_cluster_prec>::base_list[4] = {'A','C','G','T'};

#endif
