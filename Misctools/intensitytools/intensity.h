#ifndef SIMPLECALLER_INTENSITY
#define SIMPLECALLER_INTENSITY

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class ReadIntensity {
public:

  typedef int Base_Type;
  static const int base_a = 0;
  static const int base_t = 1;
  static const int base_g = 2;
  static const int base_c = 3;
  const static int base_count = 4;
  
  double bases[base_count];
  
  ReadIntensity() {
  }
  
  ReadIntensity(double a,double c,double g,double t) {
    bases[base_a] = a;
    bases[base_t] = t;
    bases[base_g] = g;
    bases[base_c] = c;
  }

  Base_Type max_base() const {

    double maxval = bases[0];
    Base_Type maxpos = 0;
    
    for(Base_Type n=0;n<base_count;n++) {
      if(bases[n] >= maxval) {maxval = bases[n]; maxpos=n;}
    }
    
    return maxpos;
  }

  double max_intensity() const {
    return bases[max_base()];
  }

  Base_Type sub_max_base() const {
    Base_Type max_position = max_base();
   

    // The next base, which isn't this one (%4 loops to start)
    double maxval = bases[(max_position+1)%4];
    Base_Type maxpos = base_count;
    
    for(Base_Type n=0;n<base_count;n++) {
      if((bases[n] >= maxval) && (n != max_position)) {maxval = bases[n]; maxpos=n;}
    }

    return maxpos;
  }

  double average() const {
    double total=0;
    for(int n=0;n<base_count;n++) {
      total += bases[n];
    }
  
    return total/base_count;
  }

  double purity() const {
    double max_position     = getbase(max_base());
    double sub_max_position = getbase(sub_max_base());
    
    if(max_position < 0) return 0;
    if(sub_max_position < 0) sub_max_position = 0;

    double sum = max_position + sub_max_position;

    if(sum<0) return 0;
    return max_position/sum;    
  }

  /// No bounds checking!
  inline double getbase(Base_Type position) const {
    return bases[position];
  }

  /// No bounds checking!
  inline void setbase(Base_Type position,double value) {
    bases[position] = value;
  }

  ReadIntensity & operator<<(std::istream &in) {
    in >> bases[base_a];
    in >> bases[base_c];
    in >> bases[base_g];
    in >> bases[base_t];

    return *this;
  }

  bool offedge() {
    if((bases[base_a] < 0.1) &&
       (bases[base_c] < 0.1) &&
       (bases[base_t] < 0.1) &&
       (bases[base_g] < 0.1)) return true; else return false;
  }
  
  bool anyoffedge() {
    if(((bases[base_a] > -0.001) && (bases[base_a] < 0.001)) ||
       ((bases[base_c] > -0.001) && (bases[base_c] < 0.001)) ||
       ((bases[base_t] > -0.001) && (bases[base_t] < 0.001)) ||
       ((bases[base_g] > -0.001) && (bases[base_g] < 0.001)) ) return true; else return false;
  }

};



class IntensitySequence {
public:

  IntensitySequence() {
  }

  IntensitySequence(int x_in,int y_in,int lane_in,int tile_in) : x(x_in),
                                                                 y(y_in),
								 lane(lane_in),
								 tile(tile_in) {
  }

  // Return minimum Purity of bases between start_base and end_base.
  // No bounds checking.
  double min_purity(int start_base,int end_base) {
    
    double min_purity = intensities[0].purity();
    for(int i=start_base;i <= end_base;i++) {
      if(intensities[i].purity() < min_purity) min_purity = intensities[i].purity();
    }
    return min_purity;
  }

  double average_purity() {

    double purity_sum=0;
    for(intensities_type::const_iterator i = intensities.begin();i != intensities.end();i++) {
      purity_sum += (*i).purity();
    }
    return purity_sum/intensities.size();
  }
 
  typedef vector<ReadIntensity> intensities_type;

  intensities_type intensities;
  int x;
  int y;

  int lane;
  int tile;

  bool valid;
};

#endif
