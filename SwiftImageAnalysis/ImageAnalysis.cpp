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

#include <iomanip>
#include "EuclideanDistanceMap.h"
#include "Watershed.h"
#include "Invert.h"
#include "NWThreshold.h"
#include "SwiftImageObject.h"
#include "Segmentation.h"
#include "swiftimagecluster_utils.h"
#include "MorphologicalErosion.h"
#include "ChannelOffsets.h"
#include "swiftimagecluster_utils.h"
#include "SwiftImageCluster.h"
#include "SimpleThreshold.h"
#include "MeanThreshold.h"
#include "MedianThreshold.h"
#include "Memstats.h"
#include "CommandLine.h"
#include "DuplicateFilter.h"

template<class _prec,class _threshold_prec>
inline void ImageAnalysis<_prec,_threshold_prec>::initialise() {

  CommandLine *parms = CommandLine::Instance();

  params_threshold_window             = parms->get_parm_as<int>("threshold_window");
  params_threshold                    = parms->get_parm_as<_prec>("threshold");
  
  params_correlation_method           = ChannelOffsets<uint16>::fft;
  params_correlation_subimages        = parms->get_parm_as<int>("correlation_subimages");
  params_correlation_subsubimages     = parms->get_parm_as<int>("correlation_subsubimages");
  params_correlation_cc_subimage_multiplier = parms->get_parm_as<int>("correlation_cc_subimage_multiplier");
  params_correlation_reference_cycle  = parms->get_parm_as<int>("correlation_reference_cycle");
  params_correlation_aggregate_cycle  = parms->get_parm_as<int>("correlation_aggregate_cycle");
  params_correlation_median_channels  = parms->get_parm_as<bool>("correlation_median_channels");
  params_correlation_threshold_window = parms->get_parm_as<int>("correlation_threshold_window");
  params_correlation_threshold        = parms->get_parm_as<_prec>("correlation_threshold");
  params_correlation_use_bases        = parms->get_parm_as<int>("correlation_use_bases");
  
  params_background_subtraction_enabled= parms->get_parm_as<bool>("background_subtraction_enabled");
  params_background_subtraction_window= parms->get_parm_as<int>("background_subtraction_window");

  params_segment_cycles               = parms->get_parm_as<int>("segment_cycles");

  params_crop                         = parms->get_parm_as<bool>("crop");
  params_crop_start_x                 = parms->get_parm_as<int>("crop_start_x");
  params_crop_end_x                   = parms->get_parm_as<int>("crop_end_x");
  params_crop_start_y                 = parms->get_parm_as<int>("crop_start_y");
  params_crop_end_y                   = parms->get_parm_as<int>("crop_end_y");

  params_dump_images                  = parms->get_parm_as<bool>("dump_images");
  params_watershed                    = parms->get_parm_as<bool>("watershed");

  params_remove_blended               = parms->get_parm_as<bool>("remove_blended");

  params_cluster_limit                = parms->get_parm_as<int>("max_clusters");

  params_calculate_noise              = parms->get_parm_as<bool>("calculate_noise");
  params_load_cycle                   = parms->get_parm_as<int>("load_cycle");

  channel_offsets_standard = NULL;
  channel_offsets_thresholded = NULL;
}

template<class _prec,class _threshold_prec>
vector<string> ImageAnalysis<_prec,_threshold_prec>::read_image_list(string image_filelist_filename) {
  ifstream image_filelist(image_filelist_filename.c_str());
  vector<string> list;

  if(!image_filelist.is_open()) {
    err << m_tt.str() << "ERROR in ImageAnalysis: Could not open image filelist at: " << image_filelist_filename << endl;
    return list; // I should probably throw an exception or pass a return code here
  }
  
  while(!image_filelist.eof()) {
    string current_line;

    getline(image_filelist, current_line);
  
    if(!image_filelist.eof()) {
      list.push_back(current_line);
    }
  }

  return list;
}

template<class _prec,class _threshold_prec>
bool ImageAnalysis<_prec,_threshold_prec>::load_images(int load_first    ,
                                                       bool get_reference
                                                      ) {
  // Clear out old images
  for(size_t n=0;n<images.size();n++) {
    images[n].clear();
  } 

  Memstats mem ("Image_Analysis::load_images");                    ///< Memory usage sampler

  // Check if all filelists are the same size, if not we have a problem (because we have differing cycle numbers for different channels)
  for(vector<vector<string> >::const_iterator i = image_filenames.begin()+1;i != image_filenames.end();i++) {
    if((*i).size() != (*(i-1)).size()) return false;
  }
  
  // We must have at least one channel (base)
  if(image_filenames.size() < 1) return false;

  int load_cycles = image_filenames[0].size();

  if(load_cycles > load_first) load_cycles=load_first;

  // Add filenames to images
  for(int current_cycle=0;current_cycle<load_cycles;current_cycle++) {
    for(int n=0;n<base_num;n++) {
      err << m_tt.str() << "Loading image: " << image_filenames[n][current_cycle];
      SwiftImage<uint16> si(image_filenames[n][current_cycle].c_str());
      //images[n].push_back(si._precsize());
      if(params_crop) {
        int start_x=0;
        int end_x  =0;
        int start_y=0;
        int end_y  =0;

        if(params_crop_start_x == -1) start_x = 0;                 else start_x = params_crop_start_x;
        if(params_crop_end_x   == -1) end_x   = si.image_width();  else end_x   = params_crop_end_x;

        if(params_crop_start_y == -1) start_y = 0;                 else start_y = params_crop_start_y;
        if(params_crop_end_y   == -1) end_y   = si.image_height(); else end_y   = params_crop_end_y;

        cout << " cropping " << start_x << "," << end_x << " " << start_y << "," << end_y;
        si = si.crop(start_x,
                     end_x,
                     start_y,
                     end_y);
        
        si.clear_offset(); 
      }
      cout << endl;
      
      images[n].push_back(si);
    }
  }

  // Clear out loaded filenames
  for(int n=0;n<base_num;n++) {
    image_filenames[n].erase(image_filenames[n].begin(),image_filenames[n].begin()+load_cycles);
  }


  if(get_reference == true) {
    reference_images.clear();
    for(int base=0;base<base_num;base++) {
      reference_images.push_back(SwiftImage<uint16>(0,0));
      reference_images[base] = images[base][params_correlation_reference_cycle];
    }
  }

  return true; // TODO: return false on load failure
}

template<class _prec,class _threshold_prec>
vector<SwiftImageCluster<_prec> > ImageAnalysis<_prec,_threshold_prec>::segmentation_and_registration(const vector<vector<SwiftImage<> > > &images) {
    
  vector<vector<SwiftImage<_threshold_prec> > > images_thresholded(base_num,vector<SwiftImage<_threshold_prec> >(params_segment_cycles,SwiftImage<_threshold_prec>(0,0)));

  int base_num     = images.size();

  // Threshold images to identify clusters
  NWThreshold<uint16,_threshold_prec> nwt(params_threshold_window,                      // Window size +/- this value
                                          params_threshold,                             // Threshold            
                                          60000,                                        // Foreground pixels become this
                                          0,                                            // Random sampling (0=off)
                                          NWThreshold<uint16>::mask_type_square);       // Mask type
  
  #if defined(_OPENMP)
    #pragma omp parallel for
  #endif
  for(int cycle=0;cycle<params_segment_cycles;cycle++) {
    for(int base=0;base<base_num;base++) {
      // Threshold image
      nwt.process_square(images[base][cycle],images_thresholded[base][cycle]);
    }
  }
  
  // Segmentation
  Segmentation<uint16,_threshold_prec,_prec> segmenter(params_watershed);
  vector<SwiftImageCluster<_prec> > image_clusters;
  image_clusters = segmenter.process(images_thresholded,images,params_segment_cycles); //2, set from parameter
  
  // Graceful fail if too many clusters identified
  if(image_clusters.size() < 1000) {
    err << "Too few clusters were generated, analysis can not continue. (less than 1000 clusters)" << endl;
    exit(1);
  }


  if(image_clusters.size() > static_cast<unsigned int>(params_cluster_limit)) {
  
    // We've gone over the cluster limit, we can process this data so throw an exception


    err << "Number of clusters passed limit (may be specified on command line)" << endl;
    err << "This limit prevents processing from taking an unacceptably long time, or using" << endl;
    err << "an unacceptable amount of memory, the program will now exit" << endl;

    exit(1);
  }
  
  if(image_clusters.size() == 0) {

    err << "No clusters were identified! Processing will now stop." << endl;  

    exit(1);
  }

  return image_clusters;
}

template<class _prec,class _threshold_prec>
void ImageAnalysis<_prec,_threshold_prec>::correlate_and_background_subtract() {
  
  bool unified_threshold = false;
  
  vector<vector<SwiftImage<_threshold_prec> > > images_thresholded;

  if(params_background_subtraction_enabled == false) unified_threshold = false;

  if(params_background_subtraction_window == params_correlation_threshold_window) {
    unified_threshold = true;
    // Thresholding
    NWThreshold<uint16,_threshold_prec>     n1_threshold(params_threshold_window,params_threshold,true);

    for(int base=0;base<base_num;base++) {
      images_thresholded.push_back(vector<SwiftImage<_threshold_prec> >());
      for(size_t cycle=0;cycle<images[base].size();cycle++) {
        images_thresholded[base].push_back(SwiftImage<_threshold_prec>(0,0));
        
        SwiftImage<> background_image(0,0);
        cerr << "Combined Threshold/Background Subtraction for image: " << base << " " << cycle << endl;
        n1_threshold.process_square(images[base][cycle],
                                    images_thresholded[base][cycle],
                                    &background_image);

        images[base][cycle] = images[base][cycle] - background_image;
      }
    }
  }

  if((!unified_threshold) && params_background_subtraction_enabled) {
    // Background subtraction
    MorphologicalErosion<uint16> morph_open(params_background_subtraction_window,false,MorphologicalErosion<uint16>::mask_type_square);

    
    int total_cycles = images[0].size();
    
    #if defined(_OPENMP)
      #pragma omp parallel for
    #endif
    for(int cycle=0;cycle<total_cycles;cycle++) {
      for(int base=0;base<base_num;base++) {
        err << m_tt.str() << "Background subtraction cycle " << right << setw(2) << cycle+1 << " Base " << ReadIntensity<_prec>::base_name[base] << endl;
        images[base][cycle] = images[base][cycle] - morph_open.process(images[base][cycle]);
      }
    }
  }

  // Calculate offsets
  if((channel_offsets_thresholded == NULL) && (channel_offsets_standard == NULL)) {
    if(unified_threshold) {
      err << "allocating thresholder unified" << endl;
      channel_offsets_thresholded = new ChannelOffsets<_threshold_prec>(params_correlation_method,
                                                            params_correlation_subimages,
                                                            params_correlation_subsubimages,
                                                            params_correlation_cc_subimage_multiplier,
                                                            params_correlation_reference_cycle,
                                                            params_correlation_aggregate_cycle,
                                                            params_correlation_threshold_window, 
                                                            params_correlation_threshold,
                                                            params_correlation_median_channels,
                                                            unified_threshold,
                                                            params_correlation_use_bases
                                                           );

    } else {
      err << "allocating thresholder standard" << endl;
      channel_offsets_standard = new ChannelOffsets<uint16>(params_correlation_method, 
                                                   params_correlation_subimages,
                                                   params_correlation_subsubimages,
                                                   params_correlation_cc_subimage_multiplier,
                                                   params_correlation_reference_cycle,
                                                   params_correlation_aggregate_cycle,
                                                   params_correlation_threshold_window, 
                                                   params_correlation_threshold,
                                                   params_correlation_median_channels,
                                                   unified_threshold,
                                                   params_correlation_use_bases
                                                  );
    }
  }

  if(unified_threshold) {
    channel_offsets_thresholded->set_correlation_reference_cycle(params_correlation_reference_cycle);
    channel_offsets_thresholded->process(images_thresholded);
    images_thresholded.clear();
    channel_offsets_thresholded->apply_offset(images);
  } else {
    err << " standard correlation" << endl;
    channel_offsets_standard->set_correlation_reference_cycle(params_correlation_reference_cycle);
    channel_offsets_standard->process(images);
    images_thresholded.clear();
    channel_offsets_standard->apply_offset(images);
  }

  // Dump xml
  offsets_xml += "\n<offsetmaps>";
  for(int base=0;base<static_cast<int>(images.size());base++) {
    for(size_t cycle=0;cycle<images[base].size();cycle++) {
      offsets_xml += "<offsetmap base=\"" + stringify(base) + "\" cycle=\"" + stringify(cycle) + "\">\n";
      offsets_xml += images[base][cycle].dump_offset_map_xml();
      offsets_xml += "</offsetmap>";
      offsets_xml += "\n";
    }
  }
  offsets_xml += "</offsetmaps>\n";

  // Dump to screen
  for(size_t cycle=0;cycle<images[0].size();cycle++) {
    for(int n=0;n<base_num;n++) {
      err << "Cycle: " << cycle << " base: " << n << " offsetmap:" << endl;
      images[n][cycle].dump_offset_map(err);
    }
  }
}

template<class _prec,class _threshold_prec>
void ImageAnalysis<_prec,_threshold_prec>::build_clusters(vector<Cluster<_prec> >           &clusters,
                                                          vector<SwiftImageCluster<_prec> > &image_clusters
                                                         ) {
  err << m_tt.str() << "Total Clusters: " << image_clusters.size() << endl;
  // 9. Save intensities and noise estimates to cluster objects
  
  for(typename vector<SwiftImageCluster<_prec> >::iterator i=image_clusters.begin();i != image_clusters.end();i++) {
    SwiftImagePosition<> pos = (*i).get_position(); 
    
    Cluster<_prec> c;
    ClusterPosition<> p(pos.x,pos.y);
    c.set_position(p);
    
    // Intensity
    typename Cluster<_prec>::signal_vec_type ints = (*i).get_intensity_sequence(images);
    typename Cluster<_prec>::signal_vec_type &intsig = c.signal("RAW");
    intsig.insert(intsig.begin(),ints.begin(),ints.end());
   
    // Noise
    if(params_calculate_noise) { 
      typename Cluster<_prec>::noise_vec_type nses = (*i).get_noise_sequence(images);
      typename Cluster<_prec>::noise_vec_type &nsesig = c.noise("RAW");
      nsesig.insert(nsesig.begin(),nses.begin(),nses.end());
    }

    clusters.push_back(c);
  }
}

template<class _prec,class _threshold_prec>
void ImageAnalysis<_prec,_threshold_prec>::extend_clusters(vector<Cluster<_prec> >           &clusters,
                                                           vector<SwiftImageCluster<_prec> > &image_clusters
                                                          ) {
  
  int c=0;
  for(typename vector<SwiftImageCluster<_prec> >::iterator i=image_clusters.begin();i != image_clusters.end();i++) {
    
    // Intensity
    typename Cluster<_prec>::signal_vec_type ints = (*i).get_intensity_sequence(images);
    typename Cluster<_prec>::signal_vec_type &intsig = clusters[c].signal("RAW");
    intsig.insert(intsig.end(),ints.begin(),ints.end());
   
    // Noise
    if(params_calculate_noise) { 
      typename Cluster<_prec>::noise_vec_type  nses = (*i).get_noise_sequence(images);
      typename Cluster<_prec>::noise_vec_type &nsesig = clusters[c].noise("RAW");
      nsesig.insert(nsesig.end(),nses.begin(),nses.end());
    }

    c++;
  }
}

template<class _prec,class _threshold_prec>
void ImageAnalysis<_prec,_threshold_prec>::generate(vector<Cluster<_prec> > &clusters) {
  
  load_images(params_load_cycle,true);

  vector<SwiftImageCluster<_prec> > image_clusters;
  
  generate_initial(clusters,image_clusters);

  for(;images[0].size() != 0;) {
    load_images(params_load_cycle);
    
    if(images[0].size() != 0) {
      generate_additional(clusters,image_clusters);
    }
  }
 

  if(channel_offsets_standard    != NULL) delete channel_offsets_standard;
  if(channel_offsets_thresholded != NULL) delete channel_offsets_thresholded;

  channel_offsets_standard    = NULL;
  channel_offsets_thresholded = NULL;
  images.clear();
  image_clusters.clear();
}

template<class _prec,class _threshold_prec>
void ImageAnalysis<_prec,_threshold_prec>::generate_initial(vector<Cluster<_prec> > &clusters,
                                                            vector<SwiftImageCluster<_prec> > &image_clusters) {

  Memstats mem ("Image_Analysis::generate");                // Time the whole call
  Memstats mem_misc;                                        // Timer for intermediate steps, gets reused
  
  
  // Correlate and Background Subtract
  err << "Correlation" << endl;
  correlate_and_background_subtract();
  
  // Segmentation
  err << "Segmentation" << endl;
  image_clusters = segmentation_and_registration(images);
 
  // Build clusters
  err << "Building clusters" << endl;
  build_clusters(clusters,image_clusters);

}

template<class _prec,class _threshold_prec>
void ImageAnalysis<_prec,_threshold_prec>::generate_additional(vector<Cluster<_prec> >           &clusters,
                                                               vector<SwiftImageCluster<_prec> > &image_clusters) {
  
  // Correlate and Background Subtract
  for(size_t n=0;n<images.size();n++) {
    images[n].push_back(reference_images[n]);
  }
  params_correlation_reference_cycle = images[0].size()-1;

  correlate_and_background_subtract();
  
  // Delete final cycle
  for(size_t n=0;n<images.size();n++) {
    images[n].erase(images[n].begin()+images[n].size()-1,images[n].end());
  }
 
  // Extend
  extend_clusters(clusters,image_clusters);
}
