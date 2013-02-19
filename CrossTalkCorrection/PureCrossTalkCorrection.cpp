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
_prec PureCrossTalkCorrection<_prec>::get_percentile(vector<double> c,int percentile_limit) {
  sort(c.begin(),c.end());

  return c[static_cast<int>((static_cast<double>(percentile_limit)/100)*static_cast<double>(c.size()))];
}
                                                

template<class _prec>
void PureCrossTalkCorrection<_prec>::make_bins(vector<double> &x,             /// X values
                                           vector<double> &y,             /// Y values
                                           vector<double> &xo,
                                           vector<double> &yo) { /// Bins must contain less than this many items
  for(size_t n=0;n<x.size();n++) {
    if(x[n] > y[n]) {
      xo.push_back(x[n]);
      yo.push_back(y[n]);
    }
  }

  // plotxy(xo,yo,"ctalk for regression");
}                                           

template<class _prec>
bool PureCrossTalkCorrection<_prec>::regression_arm(vector<double> &x,             /// X values
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
bool PureCrossTalkCorrection<_prec>::initialise(const vector<Cluster<_prec> > &clusters,int c_cycle,string signalid) {
  // I'm going to correct G/T and C/G independently
  // I'm going to base matrix construction on the first cycle only
  A_values.clear();
  C_values.clear();
  G_values.clear();
  T_values.clear();

  ClusterFilter_FirstOffEdge<_prec>         filter_offedge(signalid);
  

  ClusterFilter_AnyZero<_prec>         filter_any_zero(signalid,0);

  int discarded=0;
  
  int n=0;

  int clusters_per_bin = 1;
  err << "Clusters: " << clusters.size() << endl;
  err << "Clusters per bin: " << clusters_per_bin << endl;
  ClusterFilter_Erode_Crosstalk<_prec> filter_erode_ac(clusters,purecrosstalk_erode_num_bins,clusters_per_bin,ReadIntensity<_prec>::base_a,ReadIntensity<_prec>::base_c,source_signalid);// 100,10
  ClusterFilter_Erode_Crosstalk<_prec> filter_erode_gt(clusters,purecrosstalk_erode_num_bins,clusters_per_bin,ReadIntensity<_prec>::base_g,ReadIntensity<_prec>::base_t,source_signalid);// 100,10
  ClusterFilter_PurityHighest<_prec> filter_purity(clusters,12,purecrosstalk_purity_highest_how_many,source_signalid);

  // 1. `Generate filtered intensity pairs
  n=0;
  for(typename vector<Cluster<_prec> >::const_iterator i = clusters.begin();i != clusters.end();i++) {
    Cluster<_prec> filtered_cluster = (*i);

    filtered_cluster.set_valid(true);

    filtered_cluster = filter_offedge     .process_create(filtered_cluster);
    filtered_cluster = filter_erode_gt    .process(filtered_cluster);
    filtered_cluster = filter_erode_ac    .process(filtered_cluster);
    filtered_cluster = filter_any_zero    .process(filtered_cluster);
    filtered_cluster = filter_purity      .process(filtered_cluster);

    if(filtered_cluster.valid) {
      
      if((filtered_cluster.const_signal(signalid)[0].max_base() == ReadIntensity<_prec>::base_a) ||
         (filtered_cluster.const_signal(signalid)[0].max_base() == ReadIntensity<_prec>::base_c)) {
        A_values.push_back(filtered_cluster.signal(signalid)[c_cycle].get_base(ReadIntensity<_prec>::base_a));
        C_values.push_back(filtered_cluster.signal(signalid)[c_cycle].get_base(ReadIntensity<_prec>::base_c));
      }  

      if((filtered_cluster.const_signal(signalid)[0].max_base() == ReadIntensity<_prec>::base_g) ||
         (filtered_cluster.const_signal(signalid)[0].max_base() == ReadIntensity<_prec>::base_t)) {
        G_values.push_back(filtered_cluster.signal(signalid)[c_cycle].get_base(ReadIntensity<_prec>::base_g));
        T_values.push_back(filtered_cluster.signal(signalid)[c_cycle].get_base(ReadIntensity<_prec>::base_t));
      }
    } else discarded++;
    n++;
  }
  
  return true;
}
  
/// Applies the current matrix to these clusters
template<class _prec>
bool PureCrossTalkCorrection<_prec>::apply_correction(vector<Cluster<_prec> > &clusters, ///< Clusters to process
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
bool PureCrossTalkCorrection<_prec>::apply_correction(vector<Cluster<_prec> > &clusters, ///< Clusters to process
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
ReadIntensity<_prec> PureCrossTalkCorrection<_prec>::apply_correction(const ReadIntensity<_prec> &r) {
  
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
void PureCrossTalkCorrection<_prec>::apply_correction_values() {
  
  _prec t;
  
  for(size_t n=0;n<A_values.size();n++) {
    _prec a_val = A_values[n];
    _prec c_val = C_values[n];
    
    apply_correction(a_val,c_val,t,t);
    
    A_values[n] = a_val;
    C_values[n] = c_val;
  }

  for(size_t n=0;n<G_values.size();n++) {
    _prec g_val = G_values[n];
    _prec t_val = T_values[n];

    apply_correction(t,t,g_val,t_val);

    G_values[n] = g_val;
    T_values[n] = t_val;
  }

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
void PureCrossTalkCorrection<_prec>::apply_correction(_prec &a_val,_prec &c_val,_prec &g_val,_prec &t_val) {
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
