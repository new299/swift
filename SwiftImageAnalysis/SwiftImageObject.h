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

#ifndef SWIFTIMAGEANALYSIS_SWIFTIMAGEOBJECT
#define SWIFTIMAGEANALYSIS_SWIFTIMAGEOBJECT

#include <cstddef>
#include <iostream>
#include <vector>
#include "RLERun.h"
#include "SwiftImage.h"

using namespace std;

/// A SwiftImageObject represents an image feature, basically a connected set of pixels on an image
/// Pixels are run length encoded as RLERuns
template <class _prec=double>
class SwiftImageObject {
public:

  SwiftImageObject() : real(true),
                       intensity(0) {
  }


  /// Gets the intensity of this object, currently this is just the maximum intensity
  /// There maybe better ways of calculating this.
  template<class _iprec>
  _prec get_intensity(const SwiftImage<_iprec> &image,bool &onimage) const {
    
    _prec featuremax=0;
    bool first=true;

    for(unsigned int n=0;n < pixels.size();n++) {
      bool this_onimage=false;
      _prec runmax = pixels[n].max_pixel(image,this_onimage);
      if(this_onimage) onimage=true;
      if((runmax > featuremax) || first) {
        featuremax = runmax;
        first=false;
      }
    }

    return featuremax;
  }
  
  void apply_offset_map(const vector<vector<SwiftImagePosition<> > > &offset_map,
                        int image_width,
                        int image_height) {

    for(size_t n=0;n<pixels.size();n++) {
      if((pixels[n].pos.x >= 0) && (pixels[n].pos.y >= 0)) {
      
        int subimage_x = floor(static_cast<double>(pixels[n].pos.x)/(static_cast<double>(image_width) /static_cast<double>(offset_map.size())));
        int subimage_y = floor(static_cast<double>(pixels[n].pos.y)/(static_cast<double>(image_height)/static_cast<double>(offset_map[0].size())));
        
        if(subimage_x == offset_map.size()) subimage_x = offset_map.size()-1;
        if(subimage_y == offset_map[0].size()) subimage_y = offset_map[0].size()-1;
        
        // opposite shift
        pixels[n].pos.x = pixels[n].pos.x - offset_map[subimage_x][subimage_y].x;
        pixels[n].pos.y = pixels[n].pos.y - offset_map[subimage_x][subimage_y].y;
      } else {
        pixels[n].pos.x = pixels[n].pos.x - offset_map[0][0].x; //TODO: better than this if only 1 is 0.
        pixels[n].pos.y = pixels[n].pos.y - offset_map[0][0].y; 
      }
    }
  }

  void apply_offset_slope(const SwiftImagePosition<> &offset,const SwiftImagePosition<double> &slope) {
    for(size_t n=0;n<pixels.size();n++) {
      //pixels[n].pos.x = pixels[n].pos.x + offset.x;
      //pixels[n].pos.y = pixels[n].pos.y + offset.y;
      pixels[n].pos.x = pixels[n].pos.x + static_cast<int>(round(slope.x*static_cast<double>(pixels[n].pos.x))) + offset.x;
      pixels[n].pos.y = pixels[n].pos.y + static_cast<int>(round(slope.y*static_cast<double>(pixels[n].pos.y))) + offset.y;
    }
  }

  void dump(ostream &out) const {
    out << "pixels:" << endl;
    for(unsigned int r=0;r<pixels.size();r++) {
      out << r << " " << pixels[r].pos.x << "," << pixels[r].pos.y << "," << pixels[r].length << endl;
    }
    out << "intensity: " << intensity << endl;
    if(real) out << "real: true" << endl;
        else out << "real: false" << endl;
  }


  void set_image(SwiftImage<int> &image,int value) {
    for(size_t r=0;r<pixels.size();r++) {
      for(int x=pixels[r].pos.x;x<(pixels[r].pos.x+pixels[r].length);x++) {
        image(x,pixels[r].pos.y) = value;
      }
    }
  }
  
  vector<int> find_in_image(SwiftImage<int> &image) {
    vector<int> overlaps;

    for(size_t r=0;r<pixels.size();r++) {
      for(int x=pixels[r].pos.x;x<(pixels[r].pos.x+pixels[r].length);x++) {
        int pixel = image(x,pixels[r].pos.y);
        if(pixel != -1) overlaps.push_back(pixel);
      }
    }
    
    return overlaps;
  }

  bool real;
  _prec intensity;
  vector<RLERun<> > pixels;
};

#endif
