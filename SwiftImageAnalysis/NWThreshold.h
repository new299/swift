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

#ifndef SWIFTIMAGEANALYSIS_NWTHRESHOLD
#define SWIFTIMAGEANALYSIS_NWTHRESHOLD

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "SwiftImage.h"
#include "SwiftWindow.h"
#include "Timetagger.h"

using namespace std;

/// This class thresholded a grayscale SwiftImage into a binary image. It thresholds based
/// on the difference between this pixel and the maximum pixel in a window. For example
/// pixels that are 0.5*maxpixel maybe kept as foreground pixels.

template <class _prec=uint16,class _threshold_prec=uint16>
class NWThreshold {
public:
  static const int mask_type_square=0;  ///< Type to indicate a square mask around the central pixel
  static const int mask_type_circle=1;  ///< Type to indicate a cicular mask around the central pixel
  static const int mask_type_random=2;  ///< Type to indicate a randomly generated mask around the central pixel

  NWThreshold(int    window_size_in     = 2,                  ///< The size of the window around each pixel to examine
              double threshold_in       = 0.5f,               ///< The fraction of the maximum above which a pixel is foreground
              _prec  foregroundpixel_in = 1,                  ///< The value to set foreground pixels to
              int    random_sample_in   = 0,                  ///< Rather than examining all pixels we can perform random sampling to speed things up, if not 0 use this many samples.
              int    mask_type          = mask_type_square,   ///< The mask type to use (see mask_type_square etc.)
              bool   get_morphopen_in   = false               ///< Return the morphologically openned image as well
             ) : window_size(window_size_in),
                 foreground_pixel(foregroundpixel_in),
                 random_sample(random_sample_in),
                 threshold(threshold_in),
                 get_morphopen(get_morphopen_in) {

    if(mask_type==mask_type_square) mask_square();
    if(mask_type==mask_type_circle) mask_circle();
    if(mask_type==mask_type_random) mask_random(10);
  
  }


  /// This version is optimised for square masks
  void process_square(const SwiftImage<_prec>         &source,
                            SwiftImage<_threshold_prec> &dest,
                            SwiftImage<_prec>         *morphopenimg_in=NULL // I'm bad for doing this.
                     ) {
    SwiftImage<_prec> &morphopenimg = *morphopenimg_in;

    SwiftWindow<_prec> swin (window_size);
    
    SwiftImage<_prec> max_img = (swin.window_square(source, SwiftWindow<_prec>::find_max));
    SwiftImage<_prec> min_img = (swin.window_square(source, SwiftWindow<_prec>::find_min));

    vector<_prec> max_vals(max_img.get_image());
    vector<_prec> min_vals(min_img.get_image());

    cout << m_tt.str() << "process_square: max/min found" << endl;
    
    // SwiftImage<_prec> dest(source.image_width(),source.image_height());
    dest = source; //TODO: beter way of constructing Images of same dimension/offset.
    dest.clear_offset();

    
    // To store morphologically openned image
    if(get_morphopen) {
      morphopenimg = source; //TODO: beter way of constructing Images of same dimension/offset.
      morphopenimg.clear_offset();
    }

    for(int x=source.min_x(); x<source.max_x(); x++) {
      for(int y=source.min_y(); y<source.max_y(); y++) {
      
        int index = source.get_index(x,y);
        _prec minval = min_vals[index];
        _prec maxval = max_vals[index];

        // > NOT >= (would cause all zero image to be foreground)
        // TODO: We have to use dest on the next line, because we cleared the offset on dest. SwiftWindow currently doesn't cope with
        //       offsets correctly.
        if((static_cast<double>(dest(x,y))-static_cast<double>(minval)) 
            > ((static_cast<double>(maxval)-static_cast<double>(minval))*threshold)) {
          dest(x,y) = foreground_pixel;
        } else {
          dest(x,y) = 0;
        }
        if(get_morphopen) {
          morphopenimg(x,y) = minval;
        }
      }
    }


    //TODO: this is a messy hack, need to fix SwiftWindow so it can cope with offsetted images
    dest.copy_offset(source);
  }

  SwiftImage<_prec> process(const SwiftImage<_prec> &source) {
  
    SwiftImage<_prec> dest(source.image_width(),source.image_height());

    for(int x=0;x<source.image_width();x++) {
      for(int y=0;y<source.image_height();y++) {
        
        _prec minval  = 0;
        _prec maxval  = 0;
        bool minfirst = true;
        bool maxfirst = true;

        if(random_sample == 0) {
          for(int cx=x-window_size;cx<=(x+window_size);cx++) {
            for(int cy=y-window_size;cy<=(y+window_size);cy++) {
              
              // Only use if set in mask
              if(mask[cx-x+window_size][cy-y+window_size]) {
                // Only use pixel if in image
                if(source.onimage(cx,cy)) {
                    
                    if((source(cx,cy) < minval) || minfirst) {
                      minval = source(cx,cy);
                      minfirst  = false;
                    }
                    
                    if((source(cx,cy) > maxval) || maxfirst) {
                      maxval = source(cx,cy);
                      maxfirst  = false;
                    }
                }
              }
            }
          }
        } else {
          for(int n=0;n<random_sample;n++) {
            int cx=x-window_size+(rand()%(window_size*2));
            int cy=y-window_size+(rand()%(window_size*2));

            if(mask[cx-x+window_size][cy-y+window_size]) {
              if(source.onimage(cx,cy)) {
                if((source(cx,cy) < minval) || minfirst) {
                  minval = source(cx,cy);
                  minfirst  = false;
                }
                
                if((source(cx,cy) > maxval) || maxfirst) {
                  maxval = source(cx,cy);
                  maxfirst  = false;
                }
              }
            }
          }
        }

        if((static_cast<double>(source(x,y))-static_cast<double>(minval)) > ((static_cast<double>(maxval)-static_cast<double>(minval))*threshold)) { dest(x,y) = foreground_pixel;}
        else { dest(x,y) = 0; }
      }
    }
    
    return dest;
  }

private:
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

  
  int    window_size;                ///< Window size around pixel (-/+ this much)
  _prec  foreground_pixel;           ///< Set foreground pixels to this
  int    random_sample;              ///< Number of random samples to use (0 if not using)
  double threshold;                  ///< Fraction of maxpixel value above which to select foreground pixels
  vector<vector<bool> > mask;        ///< Stores the generated mask
  bool   get_morphopen;              ///< Store morph open image?
  Timetagger m_tt;                   ///< timetag generating object
  
};

#endif
