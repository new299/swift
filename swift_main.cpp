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

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "ImageAnalysis.h"
#include "MockImageAnalysis.h"
#include "CrossTalkCorrection.h"
#include "PureCrossTalkCorrection.h"
#include "PhasingCorrection.h"
#include "Cluster.h"
#include "PrbBaseCaller.h"
#include "FastaReader.h"
#include "FastqWriter.h"
#include "Fast4Writer.h"
#include "FastaSequence.h"
#include "Reporting.h"
#include "SmallAlign.h"
#include "stringify.h"
#include "Memstats.h"
#include "CommandLine.h"
#include "Timetagger.h"
#include "clusterfilter_makepositive.h"
#include "clusterfilter_normalise.h"
#include "clusterfilter_normalisecycle.h"
#include "clusterfilter_normalisecalls.h"
// #include "clusterfilter_scalequality.h"
#include "clusterfilter_heuristicfix.h"
#include "clusterfilter_purity.h"
#include "clusterfilter_purity2.h"
#include "clusterfilter_purityhighest.h"
#include "clusterfilter_opticalduplicates.h"

using namespace std;

typedef float _precision;

bool process_parameters (int argc, char **argv);
void write_sequence_files(CommandLine *parms, vector<Cluster<_precision> > &clusters); 

int main (int argc, char **argv) {
  
  cerr << "Swift " << SVN_REV << endl;
  Timetagger m_tt;

  if(!process_parameters (argc, argv)) return 0;
  CommandLine *parms = CommandLine::Instance();
 
  cerr << parms->dump_settings();

  bool have_reference = true;
  string reference_filename;
  if(parms->is_set("ref")) {
    reference_filename = (parms->get_parm("ref"));
  } else {
    have_reference     = false;
  }
  
  string images_a_filename     (parms->get_parm("img-a"));
  string images_c_filename     (parms->get_parm("img-c"));
  string images_g_filename     (parms->get_parm("img-g"));
  string images_t_filename     (parms->get_parm("img-t"));
  string raw_signals_filename  (parms->get_parm("sigs"));
  string intout_filename  (parms->get_parm("intout"));
  string corrected_intout_filename  (parms->get_parm("corrected_intout"));
  string report_filename       (parms->get_parm("report"));
  string pf_fastq              (parms->get_parm("pf"));
  string nonpf_fastq           (parms->get_parm("non-pf"));
  string tiletag               (parms->get_parm("tag"));
  
  size_t params_purity_length         (parms->get_parm_as<size_t>("purity_length"));

  bool   gnuplot; 
  if(parms->get_parm("gnuplot") == string("1")) gnuplot = true;
  else                                          gnuplot = false;
  
  bool   discard_offedge; 
  if(parms->get_parm("discard_offedge") == string("1")) discard_offedge = true;
  else                                                  discard_offedge = false;

  int params_optical_duplicates_distance    (parms->get_parm_as<int>("optical_duplicates_distance"));
  int params_optical_duplicates_mismatches  (parms->get_parm_as<int>("optical_duplicates_mismatches"));
  // int params_purity_length                  (parms->get_parm_as<int>("purity_length"));
  _precision params_purity_threshold            (parms->get_parm_as<double>("purity_threshold"));

  Memstats mem_main ("swift_main");   // this will time the whole run
  Memstats mem_misc;                  // this will be started and stopped repeatedly

  vector<Cluster<_precision> > clusters;
  string runxml;
  
  runxml += parms->dump_settings_xml();

  if(!parms->is_set("intfile")) {
    // Perform image analysis

    // Image analysis
    mem_misc.start ("new image_analysis");
    ImageAnalysis<_precision> *m_image_analyser = new ImageAnalysis<_precision>(images_a_filename,
                                                                                images_c_filename,
                                                                                images_g_filename,
                                                                                images_t_filename);
 
    mem_misc.start ("analyser->generate");           // this will stop the previous timer
    m_image_analyser->generate(clusters);

    runxml += m_image_analyser->offsets_xml;

    mem_misc.start ("delete m_image_analyser");      // this will stop the previous timer
    delete m_image_analyser;
    mem_misc.stop();

  } else {
    // Load data from int file
    ifstream intfile(parms->get_parm("intfile").c_str());

    for(int n=0;!intfile.eof();) {
      string str;
      getline(intfile,str);

      if(!intfile.eof()) {
        Cluster<_precision> c;
        c.read_gapipelinestr("RAW",str);

        clusters.push_back(c);
        if(n==0) cout << c;
        n++;
      }
    } 
  } 
  
  if(parms->is_set("intout")) {
    cout << m_tt.str() << "Saving intensity data" << endl;
    ofstream intout_file(intout_filename.c_str());
    for(vector<Cluster<_precision> >::iterator i=clusters.begin();i != clusters.end();i++) {
      intout_file << (*i).dump_gapipelinestr("RAW");
      intout_file << endl;
    }
    intout_file.close();
  }

  if(discard_offedge) {
    cout << m_tt.str() << " Discarding clusters offedge, either any cluster not full on the tile in the first cycle, or as specified at commandline" << endl;
    cout << m_tt.str() << " Note, this may not work correctly when processing from intensity files" << endl;

    // Through out anything with a 0 anywhere in the first cycle.
    ClusterFilter_FirstOffEdge<_precision> first_offedge("RAW");
    first_offedge.process(clusters);
    
    ClusterFilter_OffEdge<_precision> any_offedge("RAW",0);
    any_offedge.process(clusters);
  }

  remove_invalid_clusters(clusters);
  
  cout << m_tt.str() << " Clusters after removal: " << clusters.size() << endl;


  if(gnuplot) {
    cout << m_tt.str() << "Plotting uncorrected values" << endl;
    crosstalk_plot(clusters,"RAW",0,ReadIntensity<>::base_a,ReadIntensity<>::base_t,"Uncorrected");
    crosstalk_plot(clusters,"RAW",0,ReadIntensity<>::base_a,ReadIntensity<>::base_g,"Uncorrected");
    crosstalk_plot(clusters,"RAW",0,ReadIntensity<>::base_c,ReadIntensity<>::base_a,"Uncorrected");
    crosstalk_plot(clusters,"RAW",0,ReadIntensity<>::base_g,ReadIntensity<>::base_t,"Uncorrected");
    crosstalk_plot(clusters,"RAW",0,ReadIntensity<>::base_t,ReadIntensity<>::base_c,"Uncorrected");
    crosstalk_plot(clusters,"RAW",0,ReadIntensity<>::base_g,ReadIntensity<>::base_c,"Uncorrected");
  }

  cout << m_tt.str() << "Clustered found: " << clusters.size() << endl;

  // Crosstalk correction
  CrossTalkCorrection<_precision> m_crosstalk_correction(0,
                                                         20,
                                                         parms->get_parm_as<float>("crosstalk_slope_threshold"),
                                                         parms->get_parm_as<float>("crosstalk_lowerpercentile"),
                                                         parms->get_parm_as<float>("crosstalk_upperpercentile"),
                                                         parms->get_parm_as<int>  ("crosstalk_bin_size_required"),
                                                         parms->get_parm_as<int>  ("crosstalk_bin_threshold"),
                                                         parms->get_parm_as<int>  ("crosstalk_erode_clusters_per_bin"),
                                                         parms->get_parm_as<int>  ("crosstalk_erode_num_bins"),
                                                         "RAW",
                                                         "TALK1");
  mem_misc.start ("crosstalk correction");
  m_crosstalk_correction.process(clusters);
  mem_misc.stop();
  
  mem_misc.start ("Clear RAW clusters");
  clear_cluster_signal(clusters,"RAW");
  mem_misc.stop();
  
  // Pure Crosstalk correction
  PureCrossTalkCorrection<_precision> m_pcrosstalk_correction(0,
                                                              20,
                                                              parms->get_parm_as<float>("purecrosstalk_slope_threshold"),
                                                              parms->get_parm_as<int>("purecrosstalk_erode_num_bins"),
                                                              parms->get_parm_as<int>("purecrosstalk_purity_highest_how_many"),
                                                              "TALK1",
                                                              "CTALK");
  mem_misc.start ("pure crosstalk correction");
  m_pcrosstalk_correction.process(clusters);
  mem_misc.stop();
  

  mem_misc.start ("Clear TALK1 clusters");
  clear_cluster_signal(clusters,"TALK1");
  mem_misc.stop();
  
  // Set negative values to 0
  ClusterFilter_NegativeZero<_precision> m_clusterfilter_negativezero("CTALK","CTALK");
  mem_misc.start ("negativezero");
  m_clusterfilter_negativezero.process(clusters);
  mem_misc.stop();

  // Normalise signals 
  mem_misc.start ("normalisation");
  ClusterFilter_Normalise<_precision> m_normalisation(clusters,"CTALK","NORMALISED");
  mem_misc.start ("normalisation.process");
  m_normalisation.process(clusters);

  
  if(gnuplot) {
    cout << m_tt.str() << "Plotting corrected values" << endl;
    crosstalk_plot(clusters,"CTALK",0,ReadIntensity<>::base_a,ReadIntensity<>::base_t,"Corrected");
    crosstalk_plot(clusters,"CTALK",0,ReadIntensity<>::base_a,ReadIntensity<>::base_g,"Corrected");
    crosstalk_plot(clusters,"CTALK",0,ReadIntensity<>::base_c,ReadIntensity<>::base_a,"Corrected");
    crosstalk_plot(clusters,"CTALK",0,ReadIntensity<>::base_g,ReadIntensity<>::base_t,"Corrected");
    crosstalk_plot(clusters,"CTALK",0,ReadIntensity<>::base_t,ReadIntensity<>::base_c,"Corrected");
    crosstalk_plot(clusters,"CTALK",0,ReadIntensity<>::base_g,ReadIntensity<>::base_c,"Corrected");
  
    // Wait so the user can get a look at the crosstalk plots
    cout << m_tt.str() << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
  }
  
  mem_misc.start ("clear CTALK");
  clear_cluster_signal(clusters,"CTALK");
  mem_misc.stop();

  ClusterFilter_MakePositive<_precision> clusterfilter_makepositive(clusters,"NORMALISED","POSITIVE");
  mem_misc.start ("makepositive");
  clusterfilter_makepositive.process(clusters);
  mem_misc.start ("clear NORMALISED");
  clear_cluster_signal(clusters,"NORMALISED");
  mem_misc.stop();

  string sourceid = "POSITIVE";
  string targetid = "PHASE";

  int phasing_iterations = parms->get_parm_as<int>("phasing_iterations");
  for(int n=0;n<phasing_iterations;n++) {
    
    if(n==(phasing_iterations-1)) targetid= "FINAL";
    
    PhasingCorrection<_precision> m_phasing_correction(parms->get_parm_as<float>("phasing_threshold"),
                                                       parms->get_parm_as<float>("phasing_window"),
                                                       sourceid,
                                                       targetid);
    m_phasing_correction.process(clusters);
    
    if(sourceid.compare(targetid) != 0) {
      mem_misc.start ("clear in phasing: " + sourceid);
      clear_cluster_signal(clusters,sourceid);
      mem_misc.stop();
    }
    
    sourceid=targetid;
    targetid=sourceid;
  }
 
  clear_cluster_validity(clusters,true);

  // Remove optical duplicates
  cout << m_tt.str() << "Clusters, pre-optical duplicates removed: " << clusters.size() << endl;
  ClusterFilter_OpticalDuplicates<_precision> filter_duplicates("FINAL",params_optical_duplicates_distance,params_optical_duplicates_mismatches);
  mem_misc.start ("filter_duplicates");
  filter_duplicates.process(clusters);
  mem_misc.start ("remove_invalid_clusters");
  remove_invalid_clusters(clusters);
  mem_misc.stop();
  cout << m_tt.str() << "Clusters, optical duplicates removed: " << clusters.size() << endl;
  
  if(parms->is_set("sigs")) {
    cout << m_tt.str() << "Saving run data" << endl;
    ofstream signals_file(raw_signals_filename.c_str());
    for(vector<Cluster<_precision> >::iterator i=clusters.begin();i != clusters.end();i++) {
      signals_file << (*i);
    }
    signals_file.close();
  }
  
  if(parms->is_set("corrected_intout")) {
    cout << m_tt.str() << "Saving intensity data" << endl;
    ofstream intout_file(corrected_intout_filename.c_str());
    for(vector<Cluster<_precision> >::iterator i=clusters.begin();i != clusters.end();i++) {
      intout_file << (*i).dump_gapipelinestr("FINAL");
      intout_file << endl;
    }
    intout_file.close();
  }

  // Perform purity filtering
  //ClusterFilter_PurityAverage<_precision> pf_filter(params_purity_threshold,params_purity_length,"FINAL");
  ClusterFilter_Purity2<_precision> pf_filter(params_purity_length,params_purity_threshold,"FINAL");
  // ClusterFilter_PurityHighest<_precision> pf_filter(clusters,12,90000,"FINAL");
  pf_filter.process(clusters);

  mem_misc.start ("basecalling");
  PrbBaseCaller<_precision> m_base_caller_phas("FINAL","FINAL");
  m_base_caller_phas.process(clusters,false);
 
  mem_misc.stop();


  write_sequence_files(parms,clusters);
  
  mem_misc.start ("generate stats for reports");   // constructor generates stats
  
  int pair_break=0;
  if(parms->is_set("pair_break")) pair_break = parms->get_parm_as<int>("pair_break") - 1;   // 0-based index

  Reporting<_precision> rep(clusters,
                            have_reference,
                            reference_filename,
                            parms->is_set("pair_break"),
                            pair_break,
                            parms->get_parm_as<size_t>("align_every"),
                            tiletag
                           );
  mem_misc.start ("write reports");


  rep.write_report_file(report_filename,runxml);
  rep.write_report_human(cerr);
  mem_misc.stop();
  cout << m_tt.str() << "Swift complete, report file: " << report_filename << endl;
}



void write_sequence_files(CommandLine *parms,
                          vector<Cluster<_precision> > &clusters) {  
  
  Memstats mem_misc;                  // this will be started and stopped repeatedly
  
  string tiletag               (parms->get_parm("tag"));
  
  // Write FASTQs
  if(parms->is_set("fastq")) { 
    string fastq_prefix = parms->get_parm("fastq");

    if(!parms->is_set("pair_break")) {

      // Write FASTQs
      mem_misc.start ("write fastq files");
      FastqWriter pfFastq   (fastq_prefix + ".pf");
      FastqWriter nonpfFastq(fastq_prefix + ".nonpf");

      pfFastq.open();
      nonpfFastq.open();
      
      // PF Filter
      for(vector<Cluster<_precision> >::iterator i=clusters.begin();i != clusters.end();i++) {
        if((*i).is_valid())    pfFastq.write((*i).sequence("FINAL",tiletag));
        else                nonpfFastq.write((*i).sequence("FINAL",tiletag));
      }
      pfFastq.close();
      nonpfFastq.close();
      mem_misc.stop();
      
    } else {
      
      // Write FASTQs
      mem_misc.start ("write fastq files");
      
      int pair_break = parms->get_parm_as<int>("pair_break") - 1;     // convert to 0-based index
      
      FastqWriter pfFastq1   (fastq_prefix + ".pf.end1");
      FastqWriter pfFastq2   (fastq_prefix + ".pf.end2");

      FastqWriter nonpfFastq1(fastq_prefix + ".nonpf.end1");
      FastqWriter nonpfFastq2(fastq_prefix + ".nonpf.end2");

      pfFastq1.open();
      pfFastq2.open();
      nonpfFastq1.open();
      nonpfFastq2.open();
      
      // PF Filter
      for(vector<Cluster<_precision> >::iterator i=clusters.begin();i != clusters.end();i++) {
        
        if((*i).is_valid()) {
          pfFastq1.write((*i).sequence("FINAL",tiletag + ":end1").trim(0         ,pair_break-1));
          pfFastq2.write((*i).sequence("FINAL",tiletag + ":end2").trim(pair_break,-1          ));
        } else { 
          nonpfFastq1.write((*i).sequence("FINAL",tiletag + ":end1").trim(0         ,pair_break-1));
          nonpfFastq2.write((*i).sequence("FINAL",tiletag + ":end2").trim(pair_break,-1          ));
        }
      }

      pfFastq1.close();
      pfFastq2.close();
      
      nonpfFastq1.close();
      nonpfFastq2.close();
      mem_misc.stop();
    }
  }

  // Write FAST4s
  if(parms->is_set("fast4")) { 
    string fast4_prefix   = parms->get_parm("fast4");
    bool   fast4_textmode = parms->is_set("textmode");

    if ( ! parms->is_set("pair_break")) {

      // Write FASTQs
      mem_misc.start ("write fast4 files");
      Fast4Writer<size_t,float> pfFast4   (fast4_prefix + ".pf");
      Fast4Writer<size_t,float> nonpfFast4(fast4_prefix + ".nonpf");

      if(fast4_textmode) {
        pfFast4   .set_textmode();
        nonpfFast4.set_textmode();
      }
      
      pfFast4.open();
      nonpfFast4.open();
      
      // PF Filter
      for(vector<Cluster<_precision> >::iterator i=clusters.begin();i != clusters.end();i++) {
        if((*i).is_valid())     pfFast4.write((*i).sequence("FINAL",tiletag));
        else                 nonpfFast4.write((*i).sequence("FINAL",tiletag));
      }
      pfFast4.close();
      nonpfFast4.close();
      mem_misc.stop();
      
    } else {
      
      // Write FAST4s
      mem_misc.start ("write fast4 files");
      
      int pair_break = parms->get_parm_as<int>("pair_break") - 1;     // convert to 0-based index
      
      Fast4Writer<size_t,float> pfFast41   (fast4_prefix + ".pf.end1");
      Fast4Writer<size_t,float> pfFast42   (fast4_prefix + ".pf.end2");

      Fast4Writer<size_t,float>  nonpfFast41(fast4_prefix + ".nonpf.end1");
      Fast4Writer<size_t,float>  nonpfFast42(fast4_prefix + ".nonpf.end2");
      
      if(fast4_textmode) {
        pfFast41   .set_textmode();
        nonpfFast41.set_textmode();
        pfFast42   .set_textmode();
        nonpfFast42.set_textmode();
      }

      pfFast41.open();
      pfFast42.open();
      nonpfFast41.open();
      nonpfFast42.open();
      
      for(vector<Cluster<_precision> >::iterator i=clusters.begin();i != clusters.end();i++) {
        
        if((*i).is_valid()) {
          pfFast41.write((*i).sequence("FINAL",tiletag + ":end1").trim(0         ,pair_break-1));
          pfFast42.write((*i).sequence("FINAL",tiletag + ":end2").trim(pair_break,-1          ));
        } else { 
          nonpfFast41.write((*i).sequence("FINAL",tiletag + ":end1").trim(0         ,pair_break-1));
          nonpfFast42.write((*i).sequence("FINAL",tiletag + ":end2").trim(pair_break,-1          ));
        }
      }

      pfFast41.close();
      pfFast42.close();
      
      nonpfFast41.close();
      nonpfFast42.close();
      mem_misc.stop();
    }
  }
  
}

bool process_parameters (int argc, char **argv) {
  
  CommandLine *parms = CommandLine::Instance();

  parms->add_valid_parm("config"                               ,"Configuration file");
  parms->add_valid_parm("img-a"                                ,"A images filenames list");
  parms->add_valid_parm("img-c"                                ,"C images filenames list");
  parms->add_valid_parm("img-g"                                ,"G images filenames list");
  parms->add_valid_parm("img-t"                                ,"T images filenames list");
  parms->add_valid_parm("ref"                                  ,"Reference sequence for alignment (optional)");
  parms->add_valid_parm("sigs"                                 ,"Signals file (optional)");
  parms->add_valid_parm("intout"                               ,"GAPipeline style intensity file (optional)");
  parms->add_valid_parm("corrected_intout"                     ,"GAPipeline style intensity file, corrected signal (optional)");
  parms->add_valid_parm("report"                               ,"Report file");
  parms->add_valid_parm("fastq"                                ,"fastq file prefix");
  parms->add_valid_parm("fast4"                                ,"fast4 file prefix");
  parms->add_valid_parm("textmode"                             ,"write fast4 ascii files in text mode, full probabilities not ascii encoded quality scores");
  
  parms->add_valid_parm("intfile"                              ,"Load for Solexa style intensity file, instead of performing image analysis");
  parms->add_valid_parm("tag"                                  ,"Tag to write at the top of the report, for example Run ID, lane and tile (optional)");

  parms->add_valid_parm("optical_duplicates_distance"          ,"Distance in which to search of optical duplicates", false, "0");
  parms->add_valid_parm("optical_duplicates_mismatches"        ,"Number of mismatches allowed in sequences found", false, "10");
  parms->add_valid_parm("purity_length"                        ,"Length over which purity is calculated", false, "25");        
  parms->add_valid_parm("purity_threshold"                     ,"Threshold on purity (minimum value must be more than this)", false, "0.6");
  
  parms->add_valid_parm("threshold_window"                     ,"Thresholding for segmentation, distance", false, "6");
  parms->add_valid_parm("threshold"                            ,"Thresholding for segmentation, fraction of max value", false, "0.8");
  parms->add_valid_parm("watershed"                            ,"Apply watershed segmentation to thresholded images?", false, "false");

  parms->add_valid_parm("correlation_threshold_window"         ,"Thresholding used when correlating images (window size)", false, "12");
  parms->add_valid_parm("correlation_threshold"                ,"Thresholding used when correlating images (fraction of max value)", false, "0.5");
  parms->add_valid_parm("correlation_subimages"                ,"The number of subimages to use, e.g. 2 will produce different offsets for each image quarter", false, "6"); 
  parms->add_valid_parm("correlation_subsubimages"             ,"The number of subsubimages, these are used to make the offset calculation for subimages more robust", false, "2");
  parms->add_valid_parm("correlation_cc_subimage_multiplier"   ,"A multiplier for subimages, this allows crosschannel offsets to be calculated at a higher resolution to channel offsets",false,"1");  
  parms->add_valid_parm("correlation_reference_cycle"          ,"Reference cycle for first pass of correlation", false, "5");
  parms->add_valid_parm("correlation_aggregate_cycle"          ,"Sum images to this cycle for generating cross-channel offset", false, "34");
  parms->add_valid_parm("correlation_median_channels"          ,"Take the median between channels before cross-channel offseting", false, "true");
  parms->add_valid_parm("correlation_use_bases"                ,"Use how many channel in cross-correlation channel offsets?",false,"-1");
 
  parms->add_valid_parm("crosstalk_slope_threshold"            ,"Crosstalk will be iteratively corrected until the slope of the arms is less than this value",false,"0.0001");
  parms->add_valid_parm("crosstalk_lowerpercentile"            ,"Lower Percentile used in first round crosstalk correction",false,"15");
  parms->add_valid_parm("crosstalk_upperpercentile"            ,"Upper Percentile used in first round crosstalk correction",false,"85");
  parms->add_valid_parm("crosstalk_bin_size_required"          ,"Number of elements in bins during crosstalk correction"   ,false,"200");
  parms->add_valid_parm("crosstalk_bin_threshold"              ,"How close do we need to get to crosstalk_bin_size_required?",false,"80");
  parms->add_valid_parm("crosstalk_erode_clusters_per_bin"     ,"How many cluster must be in each bin from them to be retained in crosstalk correcion?",false,"1");
  parms->add_valid_parm("crosstalk_erode_num_bins"             ,"How many bins to we create during crosstalk erosion?",false,"500");

  parms->add_valid_parm("purecrosstalk_slope_threshold"        ,"Crosstalk will be iteratively corrected until the slope of the arms is less than this value (pure correcion)",false,"0.0001");
  parms->add_valid_parm("purecrosstalk_erode_num_bins"         ,"How many cluster must be in each bin from them to be retained in crosstalk correcion? (pure correction",false,"500");
  parms->add_valid_parm("purecrosstalk_purity_highest_how_many","Take the top N pure clusters as a sample for pure crosstalk correction",false,"10000");

  parms->add_valid_parm("phasing_threshold"                    ,"At most this much threshold should be corrected for",false,"0.5");
  parms->add_valid_parm("phasing_window"                       ,"Correct the phasing in this window around each base",false,"10");

  parms->add_valid_parm("background_subtraction_window"        ,"Background subtraction window size", false, "5");

  parms->add_valid_parm("segment_cycles"                       ,"Segment images up to this cycle (identify clusters in all these cycles)", false, "2");

  parms->add_valid_parm("dump_images"                          ,"Save images to disk at various stages of processing (unused?)", false, "false");

  parms->add_valid_parm("background_subtraction_enabled"       ,"Crop the input images?", false, "true");
  parms->add_valid_parm("crop"                                 ,"Crop the input images?", false, "false");
  parms->add_valid_parm("crop_start_x"                         ,"If cropping, start x position (or -1 for 0)", false, "-1");
  parms->add_valid_parm("crop_end_x"                           ,"If cropping, end x position (or -1 for max x)", false, "-1");
  parms->add_valid_parm("crop_start_y"                         ,"If cropping, start y position (or -1 for 0)", false, "-1");
  parms->add_valid_parm("crop_end_y"                           ,"If cropping, end y position (or -1 for max y)", false, "-1");

  parms->add_valid_parm("remove_blended"                       ,"Unused", false, "false");

  parms->add_valid_parm("pair_break"                           ,"Position of second end (first end length+1)",false,"");
  parms->add_valid_parm("max_clusters"                         ,"Maximum number of clusters to generate, will fail if more than this are created",false,"1500000"); 

  parms->add_valid_parm("calculate_noise"                      ,"Calculate noise estimates",false,"false");
  parms->add_valid_parm("align_every"                          ,"Align every Nth read",false,"50");
  parms->add_valid_parm("load_cycle"                           ,"Load and process this many images at a time (not this puts a limit on reference cycle and aggregate",false,"10");
  parms->add_valid_parm("phasing_iterations"                   ,"Number of phasing iterations",false,"3");
  parms->add_valid_parm("gnuplot"                              ,"Plot crosstalk with gnuplot",false,"false");

  parms->add_valid_parm("discard_offedge"                      ,"Discards any cluster that has fallen off the edge of the imaging area, in any cycle",false,"false");

  if(argc < 2) {
    cout << parms->usage() << endl;
    cout << "Where <X Images> is a line delimited list of tile images, in cycle order." << endl;
    cout << "Note: This binary processes a single tile at a time" << endl;
    return false;
  }

  parms->process_args (argc, argv);
  if (parms->is_set("config")) {
    parms->read_config_file(parms->get_parm("config"));
  }
  parms->check();
  
  return true;
  
}
