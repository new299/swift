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
#include "MockImageAnalysis.h"
#include <vector>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

template<class _prec>
inline void MockImageAnalysis<_prec>::initialise() {
}

template<class _prec>
vector<string> MockImageAnalysis<_prec>::read_image_list(string image_filelist_filename) {
  vector<string> s;
  return s;
}

template<class _prec>
bool MockImageAnalysis<_prec>::load_images() {
  
  return true; 
}

template<class _prec>
inline const vector<Cluster<_prec> > &MockImageAnalysis<_prec>::generate() {
  return generate_1pass();
}

template<class _prec>
inline const vector<Cluster<_prec> > &MockImageAnalysis<_prec>::generate_1pass() {
 
  vector<double> x;
  vector<double> y;

  int randfactor=100;

  for(int n=10000;n>0;n--) {
    Cluster<_prec> c;
    c.set_position(ClusterPosition<int>(n,n));
    _prec intensity1 = n+rand()%randfactor;
    _prec intensity2 = n/3+rand()%randfactor;
    c.raw_signal().push_back(ReadIntensity<_prec>(intensity1,intensity2,0,0));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    c.raw_signal().push_back(ReadIntensity<_prec>(intensity1/2,intensity2/2,0,20000));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    clusters.push_back(c);
  }

  for(int n=0;n<10000;n++) {
    Cluster<_prec> c;
    c.set_position(ClusterPosition<int>(n,n));
    _prec intensity1 = n/2+rand()%randfactor;
    _prec intensity2 = n+rand()%randfactor;
    c.raw_signal().push_back(ReadIntensity<_prec>(intensity1,intensity2,0,0));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    c.raw_signal().push_back(ReadIntensity<_prec>(intensity1/2,intensity2/2,0,20000));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    clusters.push_back(c);
  }
  
  for(int n=10000;n>0;n--) {
    Cluster<_prec> c;
    c.set_position(ClusterPosition<int>(n,n));

    _prec intensity1 = n+rand()%randfactor;
    _prec intensity2 = n/3+rand()%randfactor;
    c.raw_signal().push_back(ReadIntensity<_prec>(0,0,intensity1,intensity2));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    c.raw_signal().push_back(ReadIntensity<_prec>(20000,0,intensity1/2,intensity2/2));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    clusters.push_back(c);
  }

  for(int n=0;n<10000;n++) {
    Cluster<_prec> c;
    c.set_position(ClusterPosition<int>(n,n));

    _prec intensity1 = n/2+rand()%randfactor;
    _prec intensity2 = n/2+rand()%randfactor;
    c.raw_signal().push_back(ReadIntensity<_prec>(0,0,intensity1,intensity2));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    c.raw_signal().push_back(ReadIntensity<_prec>(20000,0,intensity1/2,intensity2/2));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    clusters.push_back(c);
  }

  for(int n=0;n<1000;n++) {
    Cluster<_prec> c;
    c.set_position(ClusterPosition<int>(n,n));
    _prec intensity1 = rand()%10000;
    _prec intensity2 = rand()%10000;
    _prec intensity3 = rand()%10000;
    _prec intensity4 = rand()%10000;
    c.raw_signal().push_back(ReadIntensity<_prec>(intensity1,intensity2,intensity3,intensity4));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    c.raw_signal().push_back(ReadIntensity<_prec>(intensity1/2,intensity2/2,intensity3/2,intensity4/2));
    c.raw_noise() .push_back(ReadIntensity<_prec>(0,0,0,0));
    
    clusters.push_back(c);
  }

  return get_clusters();
}

template<class _prec>
inline const vector<Cluster<_prec> > &MockImageAnalysis<_prec>::generate_2pass() {
  return generate_1pass();
}

template<class _prec>
void MockImageAnalysis<_prec>::average_offsets() {
}

template<class _prec> 
void MockImageAnalysis<_prec>::save_alignment_transforms(string filename) {
}

template<class _prec>
void MockImageAnalysis<_prec>::load_default_alignment_transforms(string filename) {
}

template<class _prec>
bool MockImageAnalysis<_prec>::run_image_analysis(int pass) {

  return true; // TODO: return false on failure
}
