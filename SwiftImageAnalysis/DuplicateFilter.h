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

#ifndef SWIFTIMAGEANALYSIS_DUPLICATEFILTER_H
#define SWIFTIMAGEANALYSIS_DUPLICATEFILTER_H

#include <cstddef>
#include <iostream>
#include <iomanip>
#include <vector>
#include "RLERun.h"
#include <math.h>
#include "SwiftImage.h"
#include "RunLengthEncode.h"
#include "DSets.h"
#include "SwiftImageObject.h"
#include "SwiftImageCluster.h"
#include "swiftimagecluster_utils.h"
#include "Timetagger.h"


using namespace std;


template <class _prec=uint16>
class DuplicateFilter {
public:
  
  DuplicateFilter(int      image_width_in,
                  int      image_height_in,
                  ostream &err_in=std::cerr) : image_width(image_width_in),
                                               image_height(image_height_in),
                                               err(err_in) {
  }

  void process(const vector<vector<SwiftImage<_prec> > > &images,
                     vector<SwiftImageCluster<> >        &clusters) {
    
    SwiftImage<int> lookup(image_width+200,image_height+200,-1);
    lookup.apply_offset(SwiftImagePosition<>(100,100));

    for(size_t i=0;i<clusters.size();i++) {
      for(size_t r=0;r < clusters[i].reference_position.pixels.size();r++) {
        int x_start = clusters[i].reference_position.pixels[r].pos.x;
        int y       = clusters[i].reference_position.pixels[r].pos.y;
        int length  = clusters[i].reference_position.pixels[r].length;

        bool skip_this=false;
        for(int x=x_start;(x<(x_start+length)) && (skip_this == false);x++) {
          if(lookup(x,y) != -1) {
          //  if(clusters[i].similarity(clusters[lookup(x,y)],images) > 30) {
              int other = lookup(x,y);
              if(clusters[i].purity(images) > clusters[lookup(x,y)].purity(images)) {
                clusters[other].set_image(lookup,-1);
                clusters[i]    .set_image(lookup,i);
                clusters[other].set_invalid();
              } else {
                clusters[i]    .set_image(lookup,-1);
                clusters[other].set_image(lookup,other);
                clusters[i]    .set_invalid();
                skip_this=true;
              }
          //  }
          } else {
            lookup(x,y) = i;
          }
        }
      }
    }
  }
  
private:

  int          image_width;       ///< Image width (used in constructing lookup images)
  int          image_height;      ///< Image height (used in constructing lookup images)
  ostream      &err;              ///< Write errors here

  Timetagger m_tt;

};  // end of class definition

#endif
