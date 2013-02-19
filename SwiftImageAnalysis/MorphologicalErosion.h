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

#ifndef SWIFTIMAGEANALYSIS_MORPHOLOGICALEROSION
#define SWIFTIMAGEANALYSIS_MORPHOLOGICALEROSION

#include <iostream>
#include <vector>
#include <math.h>
#include "SwiftImage.h"
#include "SwiftWindow.h"

using namespace std;

template <class _prec=uint16>
class MorphologicalErosion {
public:
  static const int mask_type_square=0;
  static const int mask_type_circle=1;
  static const int mask_type_random=2;

  MorphologicalErosion(int window_size_in=4,
                      bool use_second_in=false,
                      bool mask_type_in=mask_type_square) 
                     : window_size(window_size_in),
                       use_second(use_second_in),
                       mask_type(mask_type_in) {
    
    if(mask_type==mask_type_square) mask_square();
    if(mask_type==mask_type_circle) mask_circle();
    if(mask_type==mask_type_random) mask_random(10);
  
  }

  bool mask_circle() {
    // Generate mask image, a circular mask
    
    mask.clear();
    mask.insert(mask.begin(),(window_size*2)+1,vector<bool>((window_size*2)+1,false));

    int radius = window_size;

    // Outline
    for(int x=0-window_size;x<=window_size;x++) {
      int y1 = static_cast<int>(sqrt((radius*radius)-(x*x)));
      int y2 = -1*y1;
      mask[x+window_size][y1+window_size]=true;
      mask[x+window_size][y2+window_size]=true;
    }
    
    //Fill it in
    bool set=false;
    for(int x=0;x<(window_size*2)+1;x++) {
      int true_count=0;
      for(int y=0;y<(window_size*2)+1;y++) {
        if(mask[x][y]==true) true_count++;
      }

      if(true_count == 2) {
        for(int y=0;y<(window_size*2)+1;y++) {
          if((mask[x][y] == true) && (set == false)) set = true;
          else if((mask[x][y] == true) && (set == true )) set = false;

          if(set == true) mask[x][y] = true;
        }
      }
      set=false;
    }
  
    return true;
  }

  bool mask_square() {
    // Generate mask image, a square mask
    
    mask.clear();
    mask.insert(mask.begin(),(window_size*2)+1,vector<bool>((window_size*2)+1,true));
  
    return true;
  }

  bool mask_random(int num) {
    // Generate mask image, a random mask
    
    mask.clear();
    mask.insert(mask.begin(),(window_size*2)+1,vector<bool>((window_size*2)+1,false));

    for(int n=0;n<num;n++) {
      int x=rand()%(window_size*2);
      int y=rand()%(window_size*2);
      
      mask[x][y] = true;
    }

    return true;
  }

  SwiftImage<_prec> process(const SwiftImage<_prec> &source) {
    
    SwiftImage<_prec> dest = source;
    dest.clear_offset();
    if (mask_type != mask_type_square || use_second) {
      dest = (process_generic(source));    // do it the hard way
    } else {
      
      SwiftWindow<> swin (window_size);
      dest = (swin.window_square(source, SwiftWindow<>::find_min)); 
    }

    dest.copy_offset(source);

    return dest;
  }
    
  SwiftImage<_prec> process_generic (const SwiftImage<_prec> &source) {
  
    SwiftImage<_prec> dest(source.image_width(),source.image_height());

    // For each pixel, examine a window_size by window_size window around the pixel,
    // Set this pixel to the minimum value found.

    for(int x=0;x<source.image_width();x++) {
      for(int y=0;y<source.image_height();y++) {
       
        _prec oldmin = source(x,y);
        _prec minval = source(x,y);

        for(int cx=x-window_size;cx<=(x+window_size);cx++) {
          for(int cy=y-window_size;cy<=(y+window_size);cy++) {
            
            // Only use if set in mask
            if(mask[cx-x+window_size][cy-y+window_size]) {
              // Only use pixel if in image
              if(((cx >= 0) && (cx < source.image_width())) &&
                 ((cy >= 0) && (cy < source.image_height()))) {
                if(source(cx,cy) < minval) {
                  oldmin = minval;
                  minval = source(cx,cy);
                }
              }
            }
          }
        }

        // Minval has been found set in dest image
        if(!use_second) dest(x,y) = minval;
                   else dest(x,y) = oldmin;
      }
    }
    
    return dest;
  }

private:
  int  window_size;
  bool use_second;
  int  mask_type;
  vector<vector<bool> > mask;
};

#endif
