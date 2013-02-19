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

#ifndef SWIFT_CROSSTALKCORRECTION_H
#define SWIFT_CROSSTALKCORRECTION_H

#include "Cluster.h"
#include <iostream>
#include <vector>
#include "Timetagger.h"
#include "plotfunctions.h"

using namespace std;

template<class _prec=double>
class CrossTalkCorrection {
public:

  CrossTalkCorrection(int correction_cycle_in,                            ///< Build crosstalk correction based on this cycle
                      int iteration_threshold_in,                   ///< Keep trying at most this number of times.
                      _prec slope_threshold_in ,                   ///< Iteration stops when slopes are less than this
                      _prec crosstalk_lowerpercentile_in,
                      _prec crosstalk_upperpercentile_in,
                      int   crosstalk_bin_size_required_in,
                      int   crosstalk_bin_threshold_in,
                      int   crosstalk_erode_clusters_per_bin_in,
                      int   crosstalk_erode_num_bins_in,
                      string source_signalid_in,                 ///< The signal ID to get intensities from
                      string target_signalid_in, ///< The signal ID to write intensties to
                      ostream &err_in=std::cerr                        ///< Write debug/errors here
                      ): correction_cycle(correction_cycle_in),
                         iteration_threshold(iteration_threshold_in),
                         slope_threshold(slope_threshold_in),
                         crosstalk_lowerpercentile(crosstalk_lowerpercentile_in),
                         crosstalk_upperpercentile(crosstalk_upperpercentile_in),
                         crosstalk_bin_size_required(crosstalk_bin_size_required_in),
                         crosstalk_bin_threshold(crosstalk_bin_threshold_in),
                         crosstalk_erode_clusters_per_bin(crosstalk_erode_clusters_per_bin_in),
                         crosstalk_erode_num_bins(crosstalk_erode_num_bins_in),
                         source_signalid(source_signalid_in),
                         target_signalid(target_signalid_in),
                         err(err_in) {
  }

  bool process(vector<Cluster<_prec> > &clusters) {
   
    ac_m=slope_threshold+1;
    ca_m=slope_threshold+1;
    gt_m=slope_threshold+1;
    tg_m=slope_threshold+1;

    string local_source_signalid = source_signalid;
    string local_target_signalid = target_signalid;
    
    initialise(clusters,correction_cycle,source_signalid);

    make_bins(A_values,C_values,A_AC_values,C_AC_values);
    make_bins(C_values,A_values,C_CA_values,A_CA_values);
    make_bins(G_values,T_values,G_GT_values,T_GT_values);
    make_bins(T_values,G_values,T_TG_values,G_TG_values);

    // We iteratively apply the correction until almost no correction was made.
    for(int n=0;(n<iteration_threshold) && 
               ((abs(ac_m) > slope_threshold) ||
                (abs(ca_m) > slope_threshold) ||
                (abs(gt_m) > slope_threshold) ||
                (abs(tg_m) > slope_threshold));n++) {
    
      make_bins(A_values,C_values,A_AC_values,C_AC_values);
      make_bins(C_values,A_values,C_CA_values,A_CA_values);
      make_bins(G_values,T_values,G_GT_values,T_GT_values);
      make_bins(T_values,G_values,T_TG_values,G_TG_values);

      err << m_tt.str() << "Number of clusters: " << clusters.size() << endl;
      // Get slopes
      regression_arm(A_AC_values,C_AC_values,ac_m,ac_c);
      regression_arm(C_CA_values,A_CA_values,ca_m,ca_c);
      regression_arm(G_GT_values,T_GT_values,gt_m,gt_c);
      regression_arm(T_TG_values,G_TG_values,tg_m,tg_c);
    
      apply_correction_values();
      apply_correction(clusters,local_source_signalid,local_target_signalid);
      
     // plotxy(A_values,C_values,"Crosstalk AC filtered corrected");
     // plotxy(G_values,T_values,"Crosstalk GT filtered corrected");

      err << m_tt.str() << "Crosstalk Correction cycle" << endl;
      err << m_tt.str() << "ac_m: " << ac_m << endl;
      err << m_tt.str() << "ca_m: " << ca_m << endl;
      err << m_tt.str() << "gt_m: " << gt_m << endl;
      err << m_tt.str() << "tg_m: " << tg_m << endl;
    
      local_source_signalid = target_signalid;
      local_target_signalid = target_signalid;
    }

    return true;
  }

  bool initialise(const vector<Cluster<_prec> > &clusters,int cycle,string signalid);

  void apply_correction_values();
  void apply_correction(_prec &a_val,_prec &c_val,_prec &g_val,_prec &t_val);
  bool apply_correction(vector<Cluster<_prec> > &clusters,string local_source_signalid,string local_target_signalid); ///< Applies the current matrix to these clusters
  bool apply_correction(vector<Cluster<_prec> > &clusters,int cycle,string local_source_signalid,string local_target_signalid); ///< Applies the current matrix to these clusters
  ReadIntensity<_prec> apply_correction(const ReadIntensity<_prec> &r);                    ///< Applies a correction matrix to a single intensity set 
  
private:


void make_bins(vector<double> &x,             /// X values
               vector<double> &y,             /// Y values
               vector<double> &xo,
               vector<double> &yo);  /// Bins must contain less than this many items
  // 1. Bin 65th -> 95th Percentile along each axis



  /// Put a regression line though the arm of a crosstalk plot and returns the slope and position it crosses the y axis.
  bool regression_arm(vector<double> &x,
                      vector<double> &y,
                      _prec         &m,
                      _prec         &c);
 
  /// Get percentile value for a vector of _precs
  _prec get_percentile(vector<double> c,int percentile_limit); 

  // Get slopes
  _prec ac_m, ac_c;                          ///< Regression values for A/C arm
  _prec ca_m, ca_c;                          ///< Regression values for C/A arm
  _prec gt_m, gt_c;                          ///< Regression values for G/T arm
  _prec tg_m, tg_c;                          ///< Regression values for T/G arm

  int correction_cycle;                      ///< Build correct from this cycle
  _prec iteration_threshold;                 ///< Try at most this many iterations
  _prec slope_threshold;                     ///< Slopes must be less than this
  
  _prec crosstalk_lowerpercentile;
  _prec crosstalk_upperpercentile;
  int   crosstalk_bin_size_required;
  int   crosstalk_bin_threshold;
  int   crosstalk_erode_clusters_per_bin;
  int   crosstalk_erode_num_bins;

  string source_signalid;
  string target_signalid;

  ostream &err;                               ///< Error output will be writen here, set to cerr in constructor default
  Timetagger m_tt;



  vector<double> A_values;
  vector<double> T_values;
  vector<double> G_values;
  vector<double> C_values;

  vector<double> A_AC_values;
  vector<double> C_AC_values;

  vector<double> A_CA_values;
  vector<double> C_CA_values;

  vector<double> G_GT_values;
  vector<double> T_GT_values;

  vector<double> G_TG_values;
  vector<double> T_TG_values;

};

// This needs to be included here because this is a template class
#include "CrossTalkCorrection.cpp"

#endif
