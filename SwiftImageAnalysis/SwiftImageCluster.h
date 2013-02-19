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

#ifndef SWIFTIMAGEANALYSIS_SWIFTIMAGECLUSTER
#define SWIFTIMAGEANALYSIS_SWIFTIMAGECLUSTER

#include <cstddef>
#include <iostream>
#include <vector>
#include "SwiftImageObject.h"
#include "ReadIntensity.h"

using namespace std;

/// This class represents a cluster, in a single cycle. Effectively this is a collection of SwiftImageObjects
/// for each base/channel there is a vector of features which are attached to this cluster, in that channel.
template<class _prec=double>
class SwiftImageCluster {

public:  
  typedef int base_type;
  static const int base_a       = 0;
  static const int base_c       = 1;
  static const int base_g       = 2;
  static const int base_t       = 3;
  static const int base_invalid = 4;
  const static int base_count = 4;

  SwiftImageCluster(const SwiftImageObject<_prec> &o) {
    reference_position = o;
    m_valid=true;
  }


  /// Return cluster position, which is currently based on the leftmost pixel of the first channel
  /// This should be fixed.
  /// TODO: Make cluster position centroidal
  SwiftImagePosition<int> get_position() {
    SwiftImagePosition<int> pos(100,100);

    pos.x = reference_position.pixels[0].pos.x;
    pos.y = reference_position.pixels[0].pos.y;
    
    return pos;
  }

  template<class _iprec>
  typename Cluster<_prec>::signal_vec_type get_intensity_sequence(const vector<vector<SwiftImage<_iprec> > > &images) const {
    
    typename Cluster<_prec>::signal_vec_type intensity_sequence;

    for(size_t cycle=0;cycle < images[0].size();cycle++) {

      ReadIntensity<_prec> r(0,0,0,0);
      for(int base=0;base<ReadIntensity<_prec>::base_count;base++) {
        
        _prec intensity=0;
        _prec maxintensity=0;
        bool first=true;
        r.bases_offedge[base] = true;
        
        bool onimage=false;
        intensity = reference_position.get_intensity(images[base][cycle],onimage);
        
        if((intensity > maxintensity) || first) {
          maxintensity = intensity;
          first=false;
        }
          
        if(onimage) {
          r.bases_offedge[base] = false;
        }
        r.set_base(base,intensity);
      }

      intensity_sequence.push_back(r);
    }

    return intensity_sequence;
  }
  
  template<class _iprec>
  int similarity(const SwiftImageCluster<_prec> &other,const vector<vector<SwiftImage<_iprec> > > &images) const {

    typename Cluster<_prec>::signal_vec_type iseq =       get_intensity_sequence(images);
    typename Cluster<_prec>::signal_vec_type mseq = other.get_intensity_sequence(images);

    int similar=0;
    for(int n=0;n<iseq.size();n++) {
      if(iseq[n].max_base() == mseq[n].max_base()) similar++;
    }
  
    return similar;
  }

  template<class _iprec>
  double purity(const vector<vector<SwiftImage<_iprec> > > &images) const {
    vector<ReadIntensity<_prec> > iseq = get_intensity_sequence(images);

    double min_p=1;
    for(int n=0;n<12;n++) {
      double p = iseq[n].purity();
      if(p < min_p) min_p = p;
    }
  
    return min_p;
  }

  void set_image(SwiftImage<int> &image,int value) {
    reference_position.set_image(image,value);
  }

  vector<int> find_in_image(SwiftImage<int> &image) {
    return reference_position.find_in_image(image);
  }

  void set_invalid() {
    m_valid=false;
  }

  void set_valid() {
    m_valid=true;
  }

  bool isvalid() {
    return m_valid;
  }
  
  /// Gets an estimation of the noise in the region containing the cluster.
  template<class _iprec>
  typename Cluster<_prec>::noise_vec_type get_noise_sequence(const vector<vector<SwiftImage<_iprec> > > &images) const {
    typename Cluster<_prec>::noise_vec_type noise_sequence;

    for(size_t cycle=0;cycle<images[0].size();cycle++) {
      ReadIntensity<_prec> r(0,0,0,0);
      // for(int base=0;base<ReadIntensity<_prec>::base_count;base++) {
        noise_sequence.push_back(r);
      // }
    }
      
    //TODO: Implement noise estimates
    return noise_sequence;
  }

  bool m_valid;
  SwiftImageObject<_prec> reference_position;
};

#endif
