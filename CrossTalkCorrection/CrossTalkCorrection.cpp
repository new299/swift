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

#include "Cluster.h"
#include <iostream>
#include <vector>
#include <gsl/gsl_fit.h>
#include "clusterfilter_negativezero.h"
#include "clusterfilter_flowcell.h"
#include "clusterfilter_erode_crosstalk.h"
#include "clusterfilter_purity.h"
#include "clusterfilter_purityhighest.h"
#include "clusterfilter_firstoffedge.h"
#include "clusterfilter_anyzero.h"
#include "clusterfilter_blobs.h"
#include "plotfunctions.h"
#include "Timetagger.h"

template<class _prec>
_prec CrossTalkCorrection<_prec>::get_percentile(vector<double> c,int percentile_limit) {
  sort(c.begin(),c.end());

  return c[static_cast<int>((static_cast<double>(percentile_limit)/100)*static_cast<double>(c.size()))];
}
                                                

template<class _prec>
void CrossTalkCorrection<_prec>::make_bins(vector<double> &x,             /// X values
                                           vector<double> &y,             /// Y values
                                           vector<double> &xo,
                                           vector<double> &yo) { /// Bins must contain less than this many items
  // 1. Bin 65th -> 95th Percentile along each axis

  // 1.1 Find 65th and 95th percentile locations

  // 1.1.3.1 Find position in count vector where we go over lower percentile limit
  // 1.1.3.2 Find position in count vector where we go over upper percentile limit
  
  _prec percentile_lower_limit_position = get_percentile(x,crosstalk_lowerpercentile); // 65
  _prec percentile_upper_limit_position = get_percentile(x,crosstalk_upperpercentile); // 95

  err << "Lower percentile Limit: " << percentile_lower_limit_position << endl;
  err << "Upper percentile Limit: " << percentile_upper_limit_position << endl;

  // 2. Take lowest value in every bin (num_bins bins)

  vector<double> bins;
  vector<double> bins_x;

  // I'm iteratively trying bin sizes, until I find one that has ten or fewer points per bin. TODO: this is inefficient, fix it.
  
  _prec bin_size_required=crosstalk_bin_size_required;
  _prec bin_threshold = crosstalk_bin_threshold;
cout << "bin_size_required: " << bin_size_required << endl;
cout << "bin_threshold    : " << bin_threshold << endl;
  _prec average_in_bins=bin_size_required+bin_threshold+1;
  bool  last_is_more  = false;
  int current_num_bins = 100;
  _prec last_num_bins = current_num_bins*2;
  int itteration_limit=20;
  int itterations=0;
  for(;(abs(average_in_bins-bin_size_required) > bin_threshold) && (itterations < itteration_limit);itterations++) {
    // 2.1 We are taking bins across the x-axis and finding the smallest y-value in those bins
    vector<bool>  bins_first(current_num_bins,true);
    vector<int>  num_in_bins(current_num_bins,0);

    bins.clear();
    bins_x.clear();
    for(int n=0;n<=current_num_bins;n++) bins  .push_back(0);
    for(int n=0;n<=current_num_bins;n++) bins_x.push_back(0);
    _prec binsize = (percentile_upper_limit_position-percentile_lower_limit_position)/current_num_bins;
    cout << "Bin size: " << binsize << endl;

    for(unsigned int n = 0;n < x.size();n++) {
      // this is very inefficient...
      for(int current_bin_number=0;current_bin_number<current_num_bins;current_bin_number++) {
        if((current_bin_number<50) || (current_bin_number>(current_num_bins-50))) {
          _prec current_bin_start = percentile_lower_limit_position+(binsize*current_bin_number);
          _prec current_bin_end   = percentile_lower_limit_position+(binsize*current_bin_number)+binsize;

         // cout << "bin start: " << current_bin_start << endl;
         // cout << "bin   end: " << current_bin_end << endl;

          if((x[n] > current_bin_start) && (x[n] <= current_bin_end)) {
            if(bins_first[current_bin_number] || (bins[current_bin_number] > y[n])) { bins  [current_bin_number] = static_cast<double>(y[n]);
                                                                                      bins_x[current_bin_number] = static_cast<double>(x[n]); }
            bins_first[current_bin_number] = false;

            num_in_bins[current_bin_number]++;
          }
        }
      }
    }

    vector<double> bins_real;
    vector<double> bins_x_real;
    int real_bins=0;
    for(size_t n=0;n<bins_first.size();n++) {
      if((bins_first[n] == false) && (((bins[n] > 1) || (bins[n] < -1)) && ((bins_x[n] > 1) || (bins_x[n] < -1)))) {
        bins_real  .push_back(bins[n]);
        bins_x_real.push_back(bins_x[n]);
        real_bins++;
      } else {
        err << "Removing empty or zero bin: " << n << endl;
      }
    }
    bins   = bins_real;
    bins_x = bins_x_real;

 
    int total=0;
    bool first=true;
    for(vector<int>::iterator i=num_in_bins.begin();i != num_in_bins.end();i++) {
      if((*i) != 0) if(((*i) < total) || first) total = (*i); // total += (*i);
      first=false;
    }
    average_in_bins = total; //static_cast<_prec>(total)/static_cast<_prec>(current_num_bins);
    err << m_tt.str() << "Current number of bins: " << current_num_bins << endl;
    err << m_tt.str() << "Average number of items per bin: " << average_in_bins << endl;

    int this_current_num_bins = current_num_bins;
    bool is_more=false;
    if(average_in_bins>bin_size_required) {current_num_bins += abs(last_num_bins - current_num_bins)/2; is_more=true;} else
                                          {current_num_bins -= abs(last_num_bins - current_num_bins)/2; is_more=false;}

    // do it again
    if(is_more == last_is_more) {
      if(average_in_bins>bin_size_required) {current_num_bins += abs(last_num_bins - current_num_bins)/2;} else
                                            {current_num_bins -= abs(last_num_bins - current_num_bins)/2;}
    }

    last_is_more=is_more;
    last_num_bins=this_current_num_bins;

    err << m_tt.str() << "Current number of bins (new): " << current_num_bins << endl;
    if(current_num_bins < 0) current_num_bins = 1;
  }
  err << m_tt.str() << "Using " << bins.size() << " bins" << endl;

  xo = bins_x;
  yo = bins;
  // plotxy(yo,xo,"Crosstalk");
}                                           

template<class _prec>
bool CrossTalkCorrection<_prec>::regression_arm(vector<double> &x,             /// X values
                                                vector<double> &y,             /// Y values
                                                _prec         &m,             /// Slope
                                                _prec         &c) {           /// Y intercept

  // regression fit
  double c0;
  double c1;
  double cov00;
  double cov01;
  double cov11;
  double sumsq;

  /* plotxy(y,x,"Crosstalk");

  // Wait so the user can get a look at the crosstalk plots
  cout << m_tt.str() << endl << "Press ENTER to continue..." << endl;

  std::cin.clear();
  std::cin.ignore(std::cin.rdbuf()->in_avail());
  std::cin.get();
*/

  gsl_fit_linear(static_cast<double *>(&x[0]),1,
                 static_cast<double *>(&y[0]),1,
                 x.size(),
                 &c0,
                 &c1,
                 &cov00,
                 &cov01,
                 &cov11,
                 &sumsq);

  m = c1;
  c = c0;

  return true;
}



/// Generate crosstalk matrix
template<class _prec>
bool CrossTalkCorrection<_prec>::initialise(const vector<Cluster<_prec> > &clusters,int c_cycle,string signalid) {
  // I'm going to correct G/T and C/G independently
  // I'm going to base matrix construction on the first cycle only
  A_values.clear();
  C_values.clear();
  G_values.clear();
  T_values.clear();

  ClusterFilter_FirstOffEdge<_prec>         filter_offedge(signalid);
  

  ClusterFilter_AnyZero<_prec>         filter_any_zero(signalid,0);

  int discarded=0;
  
  vector<Cluster<_prec> > filtered_clusters;
  vector<Cluster<_prec> > filtered_clusters2;
  int n=0;
  /*for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();(i != clusters.end()) && (n<40000);i++) {
    filtered_clusters.push_back((*i));
    n++;
  }*/

  filtered_clusters = clusters;
  
  int clusters_per_bin = crosstalk_erode_clusters_per_bin;
  err << "Clusters: " << clusters.size() << endl;
  err << "Clusters per bin: " << clusters_per_bin << endl;
  ClusterFilter_Erode_Crosstalk<_prec> filter_erode_ac(filtered_clusters,crosstalk_erode_num_bins,clusters_per_bin,ReadIntensity<_prec>::base_a,ReadIntensity<_prec>::base_c);// 100,10
  ClusterFilter_Erode_Crosstalk<_prec> filter_erode_gt(filtered_clusters,crosstalk_erode_num_bins,clusters_per_bin,ReadIntensity<_prec>::base_g,ReadIntensity<_prec>::base_t);// 100,10
  //ClusterFilter_PurityHighest<_prec> filter_purity(clusters,12,10000,source_signalid);


  // 1. `Generate filtered intensity pairs
  n=0;
  for(typename vector<Cluster<_prec> >::const_iterator i = filtered_clusters.begin();i != filtered_clusters.end();i++) {
    Cluster<_prec> filtered_cluster = (*i);

    filtered_cluster.set_valid(true);

    filtered_cluster = filter_offedge     .process_create(filtered_cluster);
    filtered_cluster = filter_erode_gt    .process(filtered_cluster);
    filtered_cluster = filter_erode_ac    .process(filtered_cluster);
    filtered_cluster = filter_any_zero    .process(filtered_cluster);
    //if(f_purity) filtered_cluster = filter_purity      .process(filtered_cluster);

    if(filtered_cluster.valid) {
      filtered_clusters2.push_back(filtered_cluster);
      
      A_values.push_back(filtered_cluster.signal(signalid)[c_cycle].get_base(ReadIntensity<_prec>::base_a));
      C_values.push_back(filtered_cluster.signal(signalid)[c_cycle].get_base(ReadIntensity<_prec>::base_c));
        
      G_values.push_back(filtered_cluster.signal(signalid)[c_cycle].get_base(ReadIntensity<_prec>::base_g));
      T_values.push_back(filtered_cluster.signal(signalid)[c_cycle].get_base(ReadIntensity<_prec>::base_t));
    } else discarded++;
    n++;
  }
  
  err << "Crosstalk, clustered not used in regression: " << discarded << endl;
  err << "Plotting uncorrected filtered values" << endl;
 // crosstalk_plot(filtered_clusters2,signalid,0,ReadIntensity<>::base_a,ReadIntensity<>::base_c,"Uncorrected filtered for regression");
 // crosstalk_plot(filtered_clusters2,signalid,0,ReadIntensity<>::base_t,ReadIntensity<>::base_g,"Uncorrected fitlered for regression");
  
  // Wait so the user can get a look at the crosstalk plots
  /* cout << m_tt.str() << endl << "Press ENTER to continue..." << endl;

  std::cin.clear();
  std::cin.ignore(std::cin.rdbuf()->in_avail());
  std::cin.get();
*/
  return true;
}
  
/// Applies the current matrix to these clusters
template<class _prec>
bool CrossTalkCorrection<_prec>::apply_correction(vector<Cluster<_prec> > &clusters, ///< Clusters to process
                                                  string local_source_signal,        ///< Source signal id
                                                  string local_target_signal         ///< Target signal id, can be the same as source
                                                 ) {

  for(typename vector<Cluster<_prec> >::iterator i = clusters.begin();i != clusters.end();i++) {
    
    typename Cluster<_prec>::signal_vec_type old_signal;
    old_signal = (*i).const_signal(local_source_signal);
    
    typename Cluster<_prec>::signal_vec_type new_signal;

    for(typename Cluster<_prec>::signal_vec_type::const_iterator j=old_signal.begin();j != old_signal.end();j++) {
      ReadIntensity<_prec> corrected = apply_correction(*j);
      new_signal.push_back(corrected);
    }

    (*i).add_signal(local_target_signal);
    (*i).signal(local_target_signal).clear();
    (*i).signal(local_target_signal) = new_signal;
   
    //TODO: inefficient
    (*i).noise(local_target_signal)  = (*i).const_noise(local_source_signal);
  }

  return true;
}

/// Applies the current matrix to these clusters, cycle based
template<class _prec>
bool CrossTalkCorrection<_prec>::apply_correction(vector<Cluster<_prec> > &clusters, ///< Clusters to process
                                                  int a_cycle,                       ///< Cycle to apply correction to
                                                  string local_source_signal,        ///< Source signal id
                                                  string local_target_signal         ///< Target signal id, can be the same as source
                                                 ) {

  for(typename vector<Cluster<_prec> >::iterator i = clusters.begin();i != clusters.end();i++) {
    
    vector<ReadIntensity<_prec> > new_signal = (*i).const_signal(local_source_signal);

    new_signal[a_cycle] = apply_correction(new_signal[a_cycle]);

    (*i).add_signal(local_target_signal);
    (*i).signal(local_target_signal).clear();
    (*i).signal(local_target_signal) = new_signal;
   
    //TODO: inefficient
    (*i).noise(local_target_signal)  = (*i).const_noise(local_source_signal);
  }

  return true;
}

template<class _prec>
ReadIntensity<_prec> CrossTalkCorrection<_prec>::apply_correction(const ReadIntensity<_prec> &r) {
  
  ReadIntensity<_prec> r1;

  _prec a_val = r.get_base(ReadIntensity<_prec>::base_a);
  _prec c_val = r.get_base(ReadIntensity<_prec>::base_c);
  _prec g_val = r.get_base(ReadIntensity<_prec>::base_g);
  _prec t_val = r.get_base(ReadIntensity<_prec>::base_t);

  apply_correction(a_val,c_val,g_val,t_val);

  r1.set_base(ReadIntensity<_prec>::base_a,a_val);
  r1.set_base(ReadIntensity<_prec>::base_c,c_val);
  r1.set_base(ReadIntensity<_prec>::base_g,g_val);
  r1.set_base(ReadIntensity<_prec>::base_t,t_val);

  return r1;
}
 
template<class _prec>
void CrossTalkCorrection<_prec>::apply_correction_values() {
  for(size_t n=0;n<A_values.size();n++) {
    _prec a_val = A_values[n];
    _prec c_val = C_values[n];
    _prec g_val = G_values[n];
    _prec t_val = T_values[n];

    apply_correction(a_val,c_val,g_val,t_val);

    A_values[n] = a_val;
    C_values[n] = c_val;
    G_values[n] = g_val;
    T_values[n] = t_val;
  }

  _prec t;
  for(size_t n=0;n<A_AC_values.size();n++) {
    _prec aval = A_AC_values[n];
    _prec cval = C_AC_values[n];
    apply_correction(aval,cval,t,t);
    A_AC_values[n] = aval;
    C_AC_values[n] = cval;
  }
  
  for(size_t n=0;n<A_CA_values.size();n++) {
    _prec aval = A_CA_values[n];
    _prec cval = C_CA_values[n];
    apply_correction(aval,cval,t,t);
    A_CA_values[n] = aval;
    C_CA_values[n] = cval;
  }

  for(size_t n=0;n<G_GT_values.size();n++) {
    _prec gval = G_GT_values[n];
    _prec tval = T_GT_values[n];
    apply_correction(t,t,gval,tval);
    G_GT_values[n] = gval;
    T_GT_values[n] = tval;
  }

  for(size_t n=0;n<G_TG_values.size();n++) {
    _prec gval = G_TG_values[n];
    _prec tval = T_TG_values[n];
    apply_correction(t,t,gval,tval);
    G_TG_values[n] = gval;
    T_TG_values[n] = tval;
  }
}

template<class _prec>
void CrossTalkCorrection<_prec>::apply_correction(_prec &a_val,_prec &c_val,_prec &g_val,_prec &t_val) {
  // Correct A/C
  _prec v = 1/(1-(ac_m*ca_m));
  _prec old_a = a_val;
  _prec old_c = c_val;
  _prec new_a = v*old_a+-1*ca_m*v*old_c;
  _prec new_c = -1*ac_m*v*old_a+v*old_c;

  a_val = new_a;
  c_val = new_c;

  // Correct G/T
  v = 1/(1-(gt_m*tg_m));
  _prec old_t = t_val;
  _prec old_g = g_val;
  _prec new_t = v*old_t+-1*gt_m*v*old_g;
  _prec new_g = -1*tg_m*v*old_t+v*old_g;
    
  t_val = new_t;
  g_val = new_g;
}
