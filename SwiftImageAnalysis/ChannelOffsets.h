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

#ifndef SWIFTIMAGEANALYSIS_CHANNELOFFSET_H
#define SWIFTIMAGEANALYSIS_CHANNELOFFSET_H

#include "SwiftFFT.h"
#include "SwiftImagePosition.h"
#include "Timetagger.h"
#include <iostream>
#include <iomanip>
#include "Memstats.h"
#include "NWThreshold.h"
#include "MedianThreshold.h"
#include "stringify.h"

class compare_deviation
{
  public:
    compare_deviation(int d) : deviation(d) {
    }

    bool operator ()(const int p1,const int p2)
    {

      int v1 = abs(deviation - p1);
      int v2 = abs(deviation - p2);
      return(v1 < v2);
    }

  int deviation;
};


template<class _prec=uint16,class _threshold_prec=uint16>
class ChannelOffsets {

public:
  
  typedef enum {
                simple,              ///< Simple matrix multiply method 
                fft,                 ///< Use fft to perform phase correlation
                validate_fft         ///< Use simple and fft, validate results
               } correlation_type; 

  ChannelOffsets(correlation_type correlation_method_in,               ///< Use FFT of naive cross-correlation
                 int              subimages,                           ///< Number of subimages to use (subimage gets different offset)
                 int              subsubimages,                        ///< Number of subsubimages (makes subimage offset calculation more robust)
                 int              cc_subimage_multiplier,              ///< multiplier for subimages for use when crosschannel offseting.
                 int              reference_cycle,                     ///< Cycle to use as reference
                 int              aggregate_cycle,                     ///< Aggregate to this cycle when performing cross-channel offseting
                 int              threshold_window,                    ///< Threshold window
                 double           threshold,                           ///< Threshold value
                 bool             median_channels,                     ///< Take median of channel offsets? (per cross-channel)
                 bool             prethresholded,                      ///< Input images are already thresholded?
                 int              use_bases=-1                         ///< Use how many bases (-1 for all)
                )         
                : crosschannel_offsets_cache(false),
                  crosschannel_offset_cache_offset1(subimages*cc_subimage_multiplier,vector<SwiftImagePosition<> >(subimages*cc_subimage_multiplier,SwiftImagePosition<>())),
                  crosschannel_offset_cache_offset2(subimages*cc_subimage_multiplier,vector<SwiftImagePosition<> >(subimages*cc_subimage_multiplier,SwiftImagePosition<>())),
                  crosschannel_offset_cache_offset3(subimages*cc_subimage_multiplier,vector<SwiftImagePosition<> >(subimages*cc_subimage_multiplier,SwiftImagePosition<>())),
                  correlation_method(correlation_method_in),
                  params_subimages(subimages),
                  params_subsubimages(subsubimages),
                  params_cc_subimage_multiplier(cc_subimage_multiplier),
                  params_reference_cycle(reference_cycle),
                  params_aggregate_cycle(aggregate_cycle),
                  params_threshold_window(threshold_window),
                  params_threshold(threshold),
                  params_median_channels(median_channels),
                  params_prethresholded(prethresholded),
                  params_use_bases(use_bases) {
    
    if(params_use_bases == -1) params_use_bases = ReadIntensity<_prec>::base_count;
  }
  
  template<class _prec1,class _prec2>
  SwiftImagePosition<> find_image_offset_fft(const SwiftImage<_prec1> &image1_in, const SwiftImage<_prec2> &image2_in) {
    
    //TODO: SwiftFFT should be able to take different types of SwiftImage
    SwiftImage<uint16> image1(0,0);
    image1 = image1_in;
    
    SwiftImage<uint16> image2(0,0);
    image2 = image2_in;

    SwiftFFT ref (
       image1.image_width(), image1.image_height(),
       image1.image_width(), image1.image_height()    // zero-fill to twice the size
    );
    
    ref.load_SwiftImage (image1);
    ref.transform_image ();
    
    // Set up a swiftFFT object to be re-used in the loop below
    SwiftFFT target (
       image2.image_width(), image2.image_height(),
       image2.image_width(), image2.image_height()    // zero-fill to twice the size
    );
    
    target.load_SwiftImage (image2, 0-image2.min_x(), 0-image2.min_y());
    target.transform_image ();
    target.cross_correlate (ref);
    SwiftImagePosition<> offset = target.get_offset();
    
    return offset;
  }
  
  template<class _prec2>
  SwiftImagePosition<> find_image_offset_fft(const SwiftFFT &ref, const SwiftImage<_prec2> &image2_in) {
    
    //TODO: SwiftFFT should be able to take different types of SwiftImage
    
    SwiftImage<uint16> image2(0,0);
    image2 = image2_in;

    // Set up a swiftFFT object to be re-used in the loop below
    SwiftFFT target (
       image2.image_width(), image2.image_height(),
       image2.image_width(), image2.image_height()    // zero-fill to twice the size
    );
    
    target.load_SwiftImage (image2, 0-image2.min_x(), 0-image2.min_y());
    target.transform_image ();
    target.cross_correlate (ref);
    SwiftImagePosition<> offset = target.get_offset();
    
    return offset;
  }
  
  vector<SwiftImage<_threshold_prec> > split_images_unsort(const SwiftImage<_threshold_prec> &image,int subimages_num,double increment) {

    vector<SwiftImage<_threshold_prec> > subimages;

    if(subimages_num == 1) {
      subimages.push_back(image);
      return subimages;
    }
    // Only processing pixels in positive range? (need to check this)

    double x_size = subimages_num;
    double y_size = subimages_num;
    int crop_size_x = static_cast<int>(image.max_x() / x_size);
    int crop_size_y = static_cast<int>(image.max_y() / y_size);

    for(double x=0;x<=(x_size-1);x=x+increment) {
      for(double y=0;y<=(y_size-1);y=y+increment) {

        int x_start   = static_cast<int>(x*static_cast<double>(crop_size_x));
        int y_start   = static_cast<int>(y*static_cast<double>(crop_size_y));
        
        subimages.push_back(image.crop(x_start,x_start+crop_size_x,y_start,y_start+crop_size_y));
        subimages[subimages.size()-1].clear_offset();

      }
    }

    return subimages;
  }

  vector<vector<SwiftImage<_threshold_prec> > > split_images(const SwiftImage<_threshold_prec> &image,int subimages_num) {

    vector<vector<SwiftImage<_threshold_prec> > > subimages(subimages_num,
                                                   vector<SwiftImage<_threshold_prec> >(subimages_num,
                                                                               SwiftImage<_threshold_prec>(0,0)));
    // Only processing pixels in positive range? (need to check this)

    if(subimages_num == 1) {
      subimages[0][0] = image;
      return subimages;
    }

    int x_size = subimages_num;
    int y_size = subimages_num;
    int crop_size_x = image.max_x() / x_size;
    int crop_size_y = image.max_y() / y_size;

    for(int x=0;x<x_size;x++) {
      for(int y=0;y<y_size;y++) {

        int x_start   = x*crop_size_x;
        int y_start   = y*crop_size_y;
        subimages[x][y] = image.crop(x_start,x_start+crop_size_x,
                                     y_start,y_start+crop_size_y);

        subimages[x][y].clear_offset(); // All Correct?
      }
    }

    return subimages;
  }


  vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > correlation_initial_split(const vector<vector<SwiftImage<_prec> > > &images) {
    vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > sub_images;
    
    // Thresholding
    NWThreshold<_prec>     n1_threshold(params_threshold_window,params_threshold);
    
    // Return vector is base,cycle,x,y.

    // Break image in to subimages
    cerr << "Split in to subimages: " << endl;
    for(int base=0;base<static_cast<int>(images.size());base++) {
      sub_images.push_back(vector<vector<vector<SwiftImage<_threshold_prec> > > >());

      for(size_t cycle=0;cycle<images[base].size();cycle++) {
        
        // Threshold the image
        SwiftImage<_threshold_prec> image(0,0);
        if(!params_prethresholded) {
          n1_threshold.process_square(images[base][cycle],image);
        } else {
          image = images[base][cycle];
          image.make_binary();
        }

        // Split it up
        sub_images[base].push_back(split_images(image,params_subimages));
      }
    }

    return sub_images;
  }

  void correlate_within_channel(vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > &sub_images) {

    for(int base=0;base<params_use_bases;base++) {
      
      cerr << "Correlating base: " << base << endl; 
      vector<vector<vector<SwiftFFT *> > > subsubimgs_ref(sub_images[base][0].size(),
             vector<vector<SwiftFFT *> >  (sub_images[base][0][0].size()));
     
      // Prepare reference FFTs
      for(size_t x=0;x<sub_images[base][0].size();x++) {
        for(size_t y=0;y<sub_images[base][0][x].size();y++) {
          vector<SwiftImage<_threshold_prec> > subsubimgs_ref_imgs = split_images_unsort(sub_images[base][params_reference_cycle][x][y] ,params_subsubimages,1);//0.5);
          for(size_t n=0;n<subsubimgs_ref_imgs.size();n++) {
              SwiftFFT *ref = new SwiftFFT(subsubimgs_ref_imgs[n].image_width(), subsubimgs_ref_imgs[n].image_height(),
                            subsubimgs_ref_imgs[n].image_width(), subsubimgs_ref_imgs[n].image_height()    // zero-fill to twice the size
                           );

              ref->load_SwiftImage (subsubimgs_ref_imgs[n]);
              ref->transform_image ();

              subsubimgs_ref[x][y].push_back(ref);
          }
        }
      }

      for(size_t cycle=0;cycle<sub_images[base].size();cycle++) {
        for(size_t x=0;x<sub_images[base][cycle].size();x++) {
          for(size_t y=0;y<sub_images[base][cycle][x].size();y++) {
            // hmm...
            sub_images[base][cycle][x][y].clear_offset();
            sub_images[base][params_reference_cycle][x][y].clear_offset();
            
            vector<int> subsub_offsets_x;
            vector<int> subsub_offsets_y;
            vector<int> subsub_score;

            // SwiftImage<_prec> m_reference = sub_images[base][params_reference_cycle][x][y];

            // cerr << "Finished building combined reference in pass1" << endl;

            vector<SwiftImage<_threshold_prec> > subsubimgs = split_images_unsort(sub_images[base][cycle][x][y],params_subsubimages,1);//0.5); // was 0.5

            for(size_t ss_n=0;ss_n<subsubimgs.size();ss_n++) {
              // Clear offsets
              subsubimgs    [ss_n].clear_offset();

              SwiftImagePosition<> offset = find_image_offset_fft(*(subsubimgs_ref[x][y][ss_n]),subsubimgs[ss_n]);
              subsubimgs[ss_n].apply_offset(offset);
              
              subsub_offsets_x.push_back(offset.x);
              subsub_offsets_y.push_back(offset.y);
            }
            
            vector<int> subsub_offsets_x_clean;
            vector<int> subsub_offsets_y_clean;
            for(size_t n=0;n<subsub_offsets_x.size();n++) { cerr << subsub_offsets_x[n] << "," << subsub_offsets_y[n] << "  "; }
            cerr << endl;

            for(size_t n=0;n<subsub_offsets_x.size();n++) {
                subsub_offsets_x_clean.push_back(subsub_offsets_x[n]);
                subsub_offsets_y_clean.push_back(subsub_offsets_y[n]);
            }

            if(subsub_offsets_x_clean.size() == 0) {
              for(size_t n=0;n<subsub_offsets_x.size();n++) {
                subsub_offsets_x_clean.push_back(subsub_offsets_x[n]);
                subsub_offsets_y_clean.push_back(subsub_offsets_y[n]);
              }
            }

            sort(subsub_offsets_x_clean.begin(),subsub_offsets_x_clean.end());
            sort(subsub_offsets_y_clean.begin(),subsub_offsets_y_clean.end());
            
            cerr << "Pass1 Using offset: " << subsub_offsets_x_clean[subsub_offsets_x_clean.size()/2] << "," << subsub_offsets_y_clean[subsub_offsets_y_clean.size()/2] << endl;
            sub_images[base][cycle][x][y].apply_offset(SwiftImagePosition<>(subsub_offsets_x_clean[subsub_offsets_x_clean.size()/2],subsub_offsets_y_clean[subsub_offsets_y_clean.size()/2])); 
          }
        }
      }
      
      for(size_t x=0;x<subsubimgs_ref.size();x++) {
        for(size_t y=0;y<subsubimgs_ref[x].size();y++) {
          for(size_t n=0;n<subsubimgs_ref[x][y].size();n++) {
              delete subsubimgs_ref[x][y][n];
          }
        }
      }
    }
  }
    
  // Median channels
  void median_channel_offsets(vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > &sub_images) {
    int base_num               = sub_images.size();
    
    for(size_t cycle=0;cycle<sub_images[0].size();cycle++) {
      for(size_t x=0;x<sub_images[0][cycle].size();x++) {
        for(size_t y=0;y<sub_images[0][cycle][x].size();y++) {
          
          vector<int> x_offsets;
          vector<int> y_offsets;

          for(int base=0;base<params_use_bases;base++) {
            x_offsets.push_back(sub_images[base][cycle][x][y].get_offset().x);
            y_offsets.push_back(sub_images[base][cycle][x][y].get_offset().y);
          }

          // The following could be used to sort using the deviation from 0,
          // at the moment it doesn't seem to help though.
          // compare_deviation comp(0);
          // sort(x_offsets.begin(),x_offsets.end(),comp);

          sort(x_offsets.begin(),x_offsets.end());
          sort(y_offsets.begin(),y_offsets.end());

          cerr << cycle << " Median offsets X: ";
          for(size_t n=0;n<x_offsets.size();n++) cerr << x_offsets[n] << " ";
          cerr << endl;

          cerr << cycle << " Median offsets Y: ";
          for(size_t n=0;n<y_offsets.size();n++) cerr << y_offsets[n] << " ";

          cerr << "******* median:(" << x_offsets[(x_offsets.size()/3)*2] << "," << y_offsets[(y_offsets.size()/3)*2] << ")" << endl;

          for(int base=0;base<base_num;base++) {
            sub_images[base][cycle][x][y].clear_offset();
            sub_images[base][cycle][x][y].apply_offset(SwiftImagePosition<>(x_offsets[(params_use_bases/2)-0.5],y_offsets[(params_use_bases/2)-0.5]));
          }
        }
      }
    }
  }

  void correlation_resplit_for_cross_channel(const vector<vector<SwiftImage<_threshold_prec> > > &images,
                                             vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > &sub_images
                                            ) {
    // This rethresholding should be optimised out.
    NWThreshold<_prec>     n1_threshold(params_threshold_window,params_threshold);
    
    int base_num               = images.size();
    
    if(params_cc_subimage_multiplier > 1) {
      // Resplit the original images based on the new size, apply the offsets calculated above to these new images.
      // resplit
      vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > sub_images2;
      // Break image in to subimages2
      cerr << "Split in to subimages2: " << endl;
      for(int base=0;base<base_num;base++) {
        sub_images2.push_back(vector<vector<vector<SwiftImage<_threshold_prec> > > >());

        for(size_t cycle=0;cycle<images[base].size();cycle++) {
          
          // Threshold the image
          SwiftImage<_threshold_prec> image(0,0);
          if(!params_prethresholded) {
            n1_threshold.process_square(images[base][cycle],image);
          } else {
            image = images[base][cycle];
            image.make_binary();
          }

          // Split it up
          sub_images2[base].push_back(split_images(image,params_subimages*params_cc_subimage_multiplier));
        }
      }

      // apply offsets
      for(size_t base=0;base<sub_images2.size();base++) {
        for(size_t cycle=0;cycle<sub_images2[base].size();cycle++) {
        
          for(size_t x=0;x<sub_images2[base][cycle].size();x++) {
            for(size_t y=0;y<sub_images2[base][cycle][x].size();y++) {
              sub_images2[base][cycle][x][y].copy_offset(sub_images[base][cycle][x/params_cc_subimage_multiplier][y/params_cc_subimage_multiplier]);
            }
          }
        }
      }
      
      sub_images = sub_images2;
    }

  }
    
  void correlate_cross_channel_use_cache(vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > &sub_images) {
    
    int base_num               = sub_images.size();
    int cycle_num              = sub_images[0].size();
    size_t subimages_cc = params_subimages * params_cc_subimage_multiplier;

    
    // This could be moved elsewhere.
    vector<vector<vector<SwiftImagePosition<> > > > crosschannel_offset_cache(base_num,vector<vector<SwiftImagePosition<> > >(subimages_cc,vector<SwiftImagePosition<> >(subimages_cc)));

    cerr << "Using cached cross channel offsets" << endl;
    for(size_t x=0;x<subimages_cc;x++) {
      for(size_t y=0;y<subimages_cc;y++) {
        crosschannel_offset_cache[1][x][y] = crosschannel_offset_cache_offset1[x][y];
        crosschannel_offset_cache[3][x][y] = crosschannel_offset_cache_offset2[x][y];
        
        crosschannel_offset_cache[2][x][y] += crosschannel_offset_cache_offset3[x][y];
        crosschannel_offset_cache[3][x][y] += crosschannel_offset_cache_offset3[x][y];
      }
    }
    
    offsetmaps.clear();
    offsetmaps = vector<vector<vector<vector<SwiftImagePosition<> > > > >(base_num,
                        vector<vector<vector<SwiftImagePosition<> > > >  (cycle_num, // num cycles
                               vector<vector<SwiftImagePosition<> > >    (subimages_cc,
                                      vector<SwiftImagePosition<> >      (subimages_cc,SwiftImagePosition<>(0,0)))));

    // Create subimage map based on these
    for(int base=0;base<base_num;base++) {
      for(size_t cycle=0;cycle<sub_images[base].size();cycle++) {
        for(size_t x=0;x<subimages_cc;x++) {
          for(size_t y=0;y<subimages_cc;y++) {
            offsetmaps[base][cycle][x][y] = sub_images[base][cycle][x/params_cc_subimage_multiplier][y/params_cc_subimage_multiplier].get_offset()
                                          + crosschannel_offset_cache[base][x][y];
          }
        }
      }
    }
  } 


  vector<vector<vector<SwiftImage<_threshold_prec> > > > build_channel_reference(const vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > &sub_images) {

    int base_num               = sub_images.size();
    // size_t subimages_cc = params_subimages * params_cc_subimage_multiplier;

    vector<vector<vector<SwiftImage<_threshold_prec> > > > 
      reference_image(base_num,vector<vector<SwiftImage<_threshold_prec> > >(sub_images[0][0].size(),
                                      vector<SwiftImage<_threshold_prec> >  (sub_images[0][0][0].size(),SwiftImage<_threshold_prec>(0,0))));

    for(int base=0;base<base_num;base++) {
      // Initialise reference images
      for(size_t x=0;x<sub_images[base][0].size();x++) {
        for(size_t y=0;y<sub_images[base][0][x].size();y++) {
          reference_image[base][x][y] = SwiftImage<_threshold_prec>(sub_images[base][0][x][y].image_width ()+200,
                                                                    sub_images[base][0][x][y].image_height()+200);
          reference_image[base][x][y].apply_offset(SwiftImagePosition<>(100,100));
        }
      }

      // For each subimage, generate a reference against this cycle
      size_t aggregate_cycle = params_aggregate_cycle;
      if(aggregate_cycle > sub_images[0].size()) {
        aggregate_cycle = sub_images[0].size()-1; 
      }

      for(size_t cycle=0;cycle<=aggregate_cycle;cycle++) {
        cerr << "CYCLE: " << cycle << " ref img max:";
        for(size_t x=0;x<sub_images[base][cycle].size();x++) {
          for(size_t y=0;y<sub_images[base][cycle][x].size();y++) {
            reference_image[base][x][y] = (reference_image[base][x][y] + sub_images[base][cycle][x][y]);
            cerr << " " << reference_image[base][x][y].max();
          }
        }
        cerr << endl;
      }

      for(size_t x=0;x<reference_image[base].size();x++) {
        for(size_t y=0;y<reference_image[base][x].size();y++) {
          
          reference_image[base][x][y] = reference_image[base][x][y].crop(sub_images[base][0][x][y].min_x(),
                                                                         sub_images[base][0][x][y].max_x(),
                                                                         sub_images[base][0][x][y].min_y(),
                                                                         sub_images[base][0][x][y].max_y()); 
          reference_image[base][x][y].clear_offset();
        }
      }
    }

    return reference_image;
  }

  void correlate_reference_images_and_apply(vector<vector<vector<SwiftImage<_threshold_prec> > > > &reference_image,
                                            vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > &sub_images) {
    // Create combined reference image
    size_t subimages_cc = params_subimages * params_cc_subimage_multiplier;
    vector<vector<SwiftImage<_prec> > > combine_reference_image(subimages_cc,vector<SwiftImage<_prec> >(subimages_cc,SwiftImage<_prec>(0,0)));
    // #if defined(_OPENMP)
    //    #pragma omp parallel for
    // #endif
    for(size_t x=0;x<reference_image[0].size();x++) {
      for(size_t y=0;y<reference_image[0][x].size();y++) {
        //SwiftImagePosition<> offset1 = reference_image[0][x][y].find_image_offset(reference_image[1][x][y],10);
        SwiftImagePosition<> offset1 = find_image_offset_fft(reference_image[0][x][y],reference_image[1][x][y]);
        reference_image[1][x][y].apply_offset(offset1);
        
        cerr << "Offset 0 1: " << offset1.x << "," << offset1.y << endl;

        //SwiftImagePosition<> offset2 = reference_image[2][x][y].find_image_offset(reference_image[3][x][y],10);
        SwiftImagePosition<> offset2 = find_image_offset_fft(reference_image[2][x][y],reference_image[3][x][y]);
        reference_image[3][x][y].apply_offset(offset2);
        
        cerr << "Offset 2 3: " << offset2.x << "," << offset2.y << endl;

        SwiftImage<_prec> join01(0,0);
        join01 = (reference_image[0][x][y] + reference_image[1][x][y]);
        SwiftImage<_prec> join23(0,0);
        join23 = (reference_image[2][x][y] + reference_image[3][x][y]);
        //SwiftImagePosition<> offset3 = join01.find_image_offset(join23,10);
        SwiftImagePosition<> offset3 = find_image_offset_fft(join01,join23);
        join23.apply_offset(offset3);
        
        cerr << "Offset 01 23: " << offset3.x << "," << offset3.y << endl;

        combine_reference_image[x][y] = (join01 + join23);
        
        for(size_t cycle=0;cycle<sub_images[1].size();cycle++) sub_images[1][cycle][x][y].apply_offset(offset1);
        for(size_t cycle=0;cycle<sub_images[3].size();cycle++) sub_images[3][cycle][x][y].apply_offset(offset2);

        for(size_t cycle=0;cycle<sub_images[2].size();cycle++) sub_images[2][cycle][x][y].apply_offset(offset3);
        for(size_t cycle=0;cycle<sub_images[3].size();cycle++) sub_images[3][cycle][x][y].apply_offset(offset3);
      
        // Cache offset

        crosschannel_offset_cache_offset1[x][y] = offset1;
        crosschannel_offset_cache_offset2[x][y] = offset2;
        crosschannel_offset_cache_offset3[x][y] = offset3;
      }
    }
  }

  void correlate_cross_channel(vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > &sub_images) {
    // size_t subimages_cc = params_subimages * params_cc_subimage_multiplier;
    int base_num               = sub_images.size();
    cerr << "Generating cross channel offsets" << endl;

    vector<vector<vector<SwiftImage<_threshold_prec> > > > reference_images = build_channel_reference(sub_images);
    
    correlate_reference_images_and_apply(reference_images,sub_images);

    // New offset map
    offsetmaps.clear();
    offsetmaps = vector<vector<vector<vector<SwiftImagePosition<> > > > >(base_num,
                        vector<vector<vector<SwiftImagePosition<> > > >  (sub_images[0].size(), // num cycles
                               vector<vector<SwiftImagePosition<> > >    (sub_images[0][0].size(),
                                      vector<SwiftImagePosition<> >      (sub_images[0][0][0].size(),SwiftImagePosition<>(0,0)))));

    // Create subimage map based on these
    for(int base=0;base<base_num;base++) {
      for(size_t cycle=0;cycle<sub_images[base].size();cycle++) {
        for(size_t x=0;x<sub_images[base][cycle].size();x++) {
          for(size_t y=0;y<sub_images[base][cycle][x].size();y++) {
            offsetmaps[base][cycle][x][y] = sub_images[base][cycle][x][y].get_offset();
          }
        }
      }
    }
  }

  void process(const vector<vector<SwiftImage<_prec> > > &images) {
   
    // OpenMP seems to like this being an int
    //int base_num               = images.size();

    // The images is divided in to subimages for correlation, this allows each subimage to get
    // a different offset. For efficiency, the following function also applies thresholding
    // if required (for efficiency).

    vector<vector<vector<vector<SwiftImage<_threshold_prec> > > > > sub_images = correlation_initial_split(images);

    // The following method correlates images within the channel and applies those offsets to the images.
    cerr << "Correlate sub_images" << endl;
    correlate_within_channel(sub_images);

    // Median channels (take median offsets, which makes things more robust).
    if(params_median_channels) median_channel_offsets(sub_images);
    
    // Cross channel correlation starts here
    if(crosschannel_offsets_cache) {
      correlate_cross_channel_use_cache(sub_images);
    } else {
      // Optionally we can resplit the images, to use a different number of subimages for cross-correlation.
      correlation_resplit_for_cross_channel(images,sub_images);
      correlate_cross_channel(sub_images);
      crosschannel_offsets_cache=true;
    }
  }


  template<class _precimg>
  void apply_offset(vector<vector<SwiftImage<_precimg> > > &images) {
    // Apply subimage map to images
    for(int base=0;base<static_cast<int>(images.size());base++) {
      for(size_t cycle=0;cycle<images[base].size();cycle++) {
        images[base][cycle].apply_offset_map(offsetmaps[base][cycle]);
      }
    }
  }

  void set_correlation_reference_cycle(size_t c) {
    params_reference_cycle = c;
  }

  vector<vector<vector<vector<SwiftImagePosition<> > > > > offsetmaps;
  
  bool crosschannel_offsets_cache;
  vector<vector<SwiftImagePosition<> > > crosschannel_offset_cache_offset1;
  vector<vector<SwiftImagePosition<> > > crosschannel_offset_cache_offset2;
  vector<vector<SwiftImagePosition<> > > crosschannel_offset_cache_offset3;
  
  
  correlation_type correlation_method;   ///< Which correlation method to use
  int    params_subimages;
  int    params_subsubimages;
  int    params_cc_subimage_multiplier;
  int    params_reference_cycle;
  int    params_aggregate_cycle;
  int    params_threshold_window;
  double params_threshold;
  bool   params_median_channels;                    ///< Take the median of the offsets on each channel (as each should have the same stage movement)
  bool   params_prethresholded;
  int    params_use_bases;
  
private:
    
  Timetagger m_tt;
  Memstats mem;
};

#endif
