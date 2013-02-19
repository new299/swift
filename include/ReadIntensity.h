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

#ifndef SWIFT_READINTENSITY
#define SWIFT_READINTENSITY

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

template<class _prec=float>
class ReadIntensity {
public:

  typedef int base_type;
  static const int base_a = 0;
  static const int base_c = 1;
  static const int base_g = 2;
  static const int base_t = 3;
  const static int base_count = 4;
  
  _prec bases[base_count];
  bool  bases_offedge[base_count];

  static const std::string base_name[]; ///< string names for bases

  ReadIntensity() {
    for(int n=0;n<base_count;n++) bases_offedge[n] = false;
  }
  
  ReadIntensity(_prec a,_prec c,_prec g,_prec t) {
    bases[base_a] = a;
    bases[base_c] = c;
    bases[base_g] = g;
    bases[base_t] = t;
    for(int n=0;n<base_count;n++) bases_offedge[n] = false;
  }
  
  ReadIntensity(_prec a,_prec c,_prec g,_prec t,bool offedge_a,bool offedge_c,bool offedge_g,bool offedge_t) {
    bases[base_a] = a;
    bases[base_c] = c;
    bases[base_g] = g;
    bases[base_t] = t;
    bases_offedge[base_a] = offedge_a;
    bases_offedge[base_c] = offedge_c;
    bases_offedge[base_g] = offedge_g;
    bases_offedge[base_t] = offedge_t;
  }

  inline base_type max_base() const {

    _prec maxval = bases[0];
    base_type maxpos = 0;
    
    for(base_type n=1;n<base_count;n++) {
      if(bases[n] >= maxval) {maxval = bases[n]; maxpos=n;}
    }
    
    return maxpos;
  }

  inline _prec maxdiff() const {
    return max_intensity()-sub_max_intensity();
  }

  inline _prec max_intensity() const {
    return bases[max_base()];
  }
  
  inline _prec sub_max_intensity() const {
    return bases[sub_max_base()];
  }

  inline base_type sub_max_base() const {
    return sub_max_base(max_base());
  }
  
  inline base_type sub_max_base(base_type max_position) const {

    // The next base, which isn't this one (%4 loops to start)
    _prec maxval = bases[(max_position+1)%4];
    base_type maxpos = (max_position+1)%4;
    
    for(base_type n=1;n<base_count-1;n++) {// added -1
      //if((bases[n] >= maxval) && (n != max_position)) {maxval = bases[n]; maxpos=n;}
      // cout << (n+max_position+1)%4 << endl;
      if(bases[(n+max_position+1)%4] >= maxval) {maxval = bases[(n+max_position+1)%4]; maxpos=(n+max_position+1)%4;}
    }

    return maxpos;
  }

  _prec average() const {
    _prec total=0;
    for(int n=0;n<base_count;n++) {
      total += bases[n];
    }
  
    return total/base_count;
  }

  _prec purity() const {
    base_type max_base_pos  = max_base();
    _prec max_intensity     = get_base(max_base_pos);
    _prec sub_max_intensity = get_base(sub_max_base(max_base_pos));
    
    if(max_intensity < 0) return 0;
    if(sub_max_intensity < 0) sub_max_intensity = 0;

    _prec sum = max_intensity + sub_max_intensity;
    
    if(sum<0) return 0;
    return max_intensity/sum;    
  }

  /// No bounds checking!
  inline _prec get_base(base_type position) const {
    return bases[position];
  }

  /// No bounds checking!
  inline void set_base(base_type position,_prec value) {
    bases[position] = value;
  }

  string as_string() const {
    string s;

    s += stringify(bases[base_a]) + " ";
    s += stringify(bases[base_c]) + " ";
    s += stringify(bases[base_g]) + " ";
    s += stringify(bases[base_t]);
  
    return s;
  }

  ReadIntensity & operator<<(std::istream &in) {
    in >> bases[base_a];
    in >> bases[base_c];
    in >> bases[base_g];
    in >> bases[base_t];

    return *this;
  }

  bool operator==(const ReadIntensity &rhs) const {
    
    if((bases[base_a] == rhs.get_base(base_a)) &&
       (bases[base_c] == rhs.get_base(base_c)) &&
       (bases[base_g] == rhs.get_base(base_g)) &&
       (bases[base_t] == rhs.get_base(base_t))) return true;
    else return false;

  }

  bool off_edge() const {
    if((bases_offedge[base_a]) &&
       (bases_offedge[base_c]) &&
       (bases_offedge[base_t]) &&
       (bases_offedge[base_g])) return true; else return false;
  }
  
  bool any_off_edge() const {
    if((bases_offedge[base_a]) ||
       (bases_offedge[base_c]) ||
       (bases_offedge[base_g]) ||
       (bases_offedge[base_t])) return true; else return false;
  }


  _prec distance(_prec pa,_prec pc,_prec pg,_prec pt,
                 _prec qa,_prec qc,_prec qg,_prec qt) const {
    return sqrt((pa-qa)*(pa-qa)+(pt-qt)*(pt-qt)+(pg-qg)*(pg-qg)+(pc-qc)*(pc-qc));
  }

  _prec get_quality() const {
    
    // probability is purity
    return (purity()-0.5)*2;
  }
};

template<class _prec>
const std::string ReadIntensity<_prec>::base_name[] = {"A", "C", "G", "T", "INVALID" };

#include "ReadIntensity.cpp"

#endif
