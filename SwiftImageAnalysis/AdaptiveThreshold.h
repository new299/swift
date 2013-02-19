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

#ifndef SWIFTIMAGEANALYSIS_ADAPTIVETHRESHOLDING
#define SWIFTIMAGEANALYSIS_ADAPTIVETHRESHOLDING

#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include <math.h>

using namespace std;

template <class _prec=uint16>
class AdaptiveThreshold {
public:
  static const int mask_type_square=0;
  static const int mask_type_circle=1;
  static const int mask_type_random=2;

  AdaptiveThreshold(int window_size_in=4,
                    double mean_offset_threshold_in=10,
                    bool mask_type=mask_type_square)
                   : window_size(window_size_in),
                     mean_offset_threshold(mean_offset_threshold_in) {

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
  
    SwiftImage<_prec> dest(source.image_width(),source.image_height());

    // For each pixel, examine a window_size by window_size window around the pixel,
    // Set this pixel to the minimum value found.

    for(int x=0;x<source.image_width();x++) {
      for(int y=0;y<source.image_height();y++) {
        
        _prec meanval = 0;
        int   numused = 0;

        for(int cx=x-window_size;cx<=(x+window_size);cx++) {
          for(int cy=y-window_size;cy<=(y+window_size);cy++) {
            
            // Only use if set in mask
            if(mask[cx-x+window_size][cy-y+window_size]) {
              // Only use pixel if in image
              if(((cx >= 0) && (cx < source.image_width())) &&
                 ((cy >= 0) && (cy < source.image_height()))) {
                  meanval += source(cx,cy);
                  numused++;
              }
            }
          }
        }

        // Minval has been found set in dest image
        meanval = meanval/numused;

        //Apply threshold
        if((source(x,y) - meanval) > meanval*mean_offset_threshold) { dest(x,y) = 65535; }
        else { dest(x,y) = 0; }
      }
    }
    
    return dest;
  }

private:
  int window_size;
  double mean_offset_threshold;
  vector<vector<bool> > mask;
};

#endif
