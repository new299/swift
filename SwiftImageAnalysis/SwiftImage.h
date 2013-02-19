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

#ifndef SWIFTIMAGEANALYSIS_SWIFTIMAGE
#define SWIFTIMAGEANALYSIS_SWIFTIMAGE

#include <iostream>
#include <vector>
#include <exception>
#include <stdexcept>
#include <iomanip>
#include <tiffio.h>
#include "SwiftImagePosition.h"
#include "Timetagger.h"
#include <sstream>
#include <math.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;

/// This class represents a SwiftImage, it currents has functions to load from and save to a tiff.
template <class _prec=uint16>
class SwiftImage {

public:
  SwiftImage(const char *filename,
             ostream &err_in=std::cerr) : image_offset_x(0),
                                          image_offset_y(0),
                                          image_slope_x(0),
                                          image_slope_y(0),
                                          offset_map(0),
                                          offset_map_enable(false),
                                          err(err_in) {
    bool rc = load(filename);
    if ( ! rc) {
        throw (std::invalid_argument("Could not load image"));
    }
    offset_map_enable=false;
  }
  
  SwiftImage(unsigned int x,                ///< X Size
             unsigned int y,                ///< Y Size
             _prec init_val=0,              ///< Initially set all pixels to this
             ostream &err_in=std::cerr      ///< Write errors/debug here
            ) : m_image_width(x),
                m_image_height(y),
                image_offset_x(0),
                image_offset_y(0),
                image_slope_x(0),
                image_slope_y(0),
                cache_max_ok(false),
                offset_map(0),
                offset_map_enable(false),
                err(err_in) {
    image.clear();

    image.insert(image.begin(),x*y,init_val);//TODO: apparently not
    offset_map_enable=false;
  }

  bool load(const char *filename);
  bool save(const char *filename) const;

  vector<_prec> &get_image() {
    return image;
  }

  const vector<_prec> &get_image() const {
    return image;
  }

  inline _prec &operator()(SwiftImagePosition<int> c) {
    cache_max_ok=false;
    bool off_edge=false;
    return const_cast<_prec &>((*this).at(c.x,c.y,off_edge)); 
  }
 
  inline  _prec &operator()(const int x, const int y) {
    cache_max_ok=false;
    bool off_edge=false;
    return const_cast<_prec &>((*this).at(x,y,off_edge));
  } 

  inline  const _prec &operator()(const int x, const int y) const {
    bool off_edge=false;
    return (*this).at(x,y,off_edge);
  }

  inline  _prec &operator()(const int x,     ///< X Position
                            const int y,     ///< Y Position
                            bool &off_edge   ///< Will return true here if the pixel is off the edge of the image TODO: Refactor, never actually set or useful (check with onimage if required), better use min_x,min_y to detect off edgeness...
                           ) {
    cache_max_ok=false;
    return const_cast<_prec &>((*this).at(x,y,off_edge));
  } 
  
  inline  const _prec &operator()(const int x, const int y,bool &off_edge) const {
    return (*this).at(x,y,off_edge);
  }

  /// This is the underlying pixel accessor all const and non-const versions call this.
  // THE AT METHOD
  inline  const _prec &at(const int x, const int y,bool &off_edge) const {
   
    int real_x=0;
    int real_y=0;

    to_real(x,y,real_x,real_y);

    if((real_x >= 0) &&
       (real_y >= 0) &&
       (real_x <  image_width()) &&
       (real_y <  image_height())) {
      return image[(real_y*m_image_width)+real_x];
    } else {
      throw out_of_range("Off edge of image");
    }
  }
  
  inline const _prec &operator()(SwiftImagePosition<int> c) const {
    bool off_edge=false;
    return (*this).at(c.x,c.y,off_edge);
  }

  SwiftImage<_prec> operator-(const SwiftImage<_prec> &rhs) const {
    SwiftImage<_prec> ret(0,0);
    
    ret = (*this);

    // Images must be of the same size
    if(rhs.image_width()  != image_width() ) return ret;
    if(rhs.image_height() != image_height()) return ret;

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        if(rhs.onimage(x,y)) {
          ret(x,y) = (*this)(x,y) - rhs(x,y);
        }
      }
    }

    return ret;
  }
  
  SwiftImage<_prec> operator/(_prec c) const {
    SwiftImage<_prec> ret(0,0);
    
    ret=(*this);

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        ret(x,y) = (*this)(x,y) / c;
      }
    }

    return ret;
  }
 

  SwiftImage<_prec> operator+(_prec c) const {
    SwiftImage<_prec> ret(0,0);
    
    ret = (*this);

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        ret(x,y) = (*this)(x,y) + c;
      }
    }

    return ret;
  }

  SwiftImage<_prec> operator+(const SwiftImage<_prec> &rhs) const {
    SwiftImage<_prec> ret = (*this);

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        if(rhs.onimage(x,y)) {
          ret(x,y) = (*this)(x,y) + rhs(x,y);
        } 
      }
    }

    return ret;
  }
  
  SwiftImage<_prec> operator&&(const SwiftImage<_prec> &rhs) const {
    SwiftImage<_prec> ret(0,0);
    ret = (*this);

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        if(rhs.onimage(x,y)) {
          if((rhs(x,y) != 0) && ((*this)(x,y) != 0)) {
            ret(x,y) = (*this)(x,y);
          } else ret(x,y) = 0;
        }
      }
    }

    return ret;
  }
  
  SwiftImage<_prec> operator/(const SwiftImage<_prec> &rhs) const {
    SwiftImage<_prec> ret(0,0);
    ret= (*this);
    
    // Images must be of the same size
    if(rhs.image_width()  != image_width() ) return ret;
    if(rhs.image_height() != image_height()) return ret;

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        if(rhs.onimage(x,y)) {
          if(rhs(x,y) != 0) {
            ret(x,y) = static_cast<_prec>(static_cast<double>((*this)(x,y)) / static_cast<double>(rhs(x,y)));
          } else ret(x,y) = 0;
        }
      }
    }

    return ret;
  }
  
  template<class _prec2>
  void copy_offset(const SwiftImage<_prec2> &other) {
    image_offset_x = other.image_offset_x;
    image_offset_y = other.image_offset_y;
    image_slope_x  = other.image_slope_x;
    image_slope_y  = other.image_slope_y;
    
    offset_map_enable = other.offset_map_enable;
    offset_map        = other.offset_map;
  }
  
  SwiftImage<_prec> &operator=(const SwiftImage<_prec> &rhs) {
    m_image_width  = rhs.m_image_width;
    m_image_height = rhs.m_image_height;

    image = rhs.image;
   
    image_offset_x = rhs.image_offset_x;
    image_offset_y = rhs.image_offset_y;
    image_slope_x  = rhs.image_slope_x;
    image_slope_y  = rhs.image_slope_y;
    
    offset_map_enable = rhs.offset_map_enable;
    offset_map        = rhs.offset_map;

    return (*this);
  }
 
  //TODO: Consolidate two version of assignment operator
  template<class _prec2>
  SwiftImage<_prec> &operator=(const SwiftImage<_prec2> &rhs) {
    m_image_width  = rhs.m_image_width;
    m_image_height = rhs.m_image_height;

    image.clear();
    for(typename vector<_prec2>::const_iterator i=rhs.image.begin();i != rhs.image.end();i++) {
      image.push_back((*i));
    }

    image_offset_x = rhs.image_offset_x;
    image_offset_y = rhs.image_offset_y;
    image_slope_x  = rhs.image_slope_x;
    image_slope_y  = rhs.image_slope_y;
    
    offset_map_enable = rhs.offset_map_enable;
    offset_map        = rhs.offset_map;

    return (*this);
  }
  
  template<class _prec2>
  SwiftImage<_prec> &copy_with_offset(const SwiftImage<_prec2> &rhs) {
    for(int x=rhs.min_x();x<rhs.max_x();x++) {
      for(int y=rhs.min_y();y<rhs.max_y();y++) {
        (*this)(x,y) = rhs(x,y);
      }
    }

    return (*this);
  }

  void make_binary(int threshold=0) {
    for(typename vector<_prec>::iterator i=image.begin();i != image.end();i++) {
      if((*i) > threshold) (*i) = 1; 
                      else (*i) = 0;
    }
  }
  
  void threshold(int thres=0) {
    for(typename vector<_prec>::iterator i=image.begin();i != image.end();i++) {
      if((*i) < thres) (*i) = 0;
    }
  }

  bool inline to_real(int x,int y,int &real_x,int &real_y) const {
    if(offset_map_enable) {
      
      if((x >= 0) && (y >= 0)) { //TODO: Why X/Y greater than 0?
        
        double dbl_subimage_x = static_cast<double>(x)
                             / (static_cast<double>(max_x())
                             /  static_cast<double>(offset_map.size()));
        double dbl_subimage_y = static_cast<double>(y)
                             / (static_cast<double>(max_y())
                             /  static_cast<double>(offset_map[0].size()));
        
        size_t subimage_x = static_cast<size_t>(floor(dbl_subimage_x));
        size_t subimage_y = static_cast<size_t>(floor(dbl_subimage_y));
        
        if(subimage_x == offset_map   .size()) subimage_x=offset_map   .size()-1;
        if(subimage_y == offset_map[0].size()) subimage_y=offset_map[0].size()-1;
        
        real_x = image_offset_x+x+offset_map[subimage_x][subimage_y].x;
        real_y = image_offset_y+y+offset_map[subimage_x][subimage_y].y;
      
      } else {
        real_x = image_offset_x+x+offset_map[0][0].x;  //TODO: do better than this if only 1 is 0.
        real_y = image_offset_y+y+offset_map[0][0].y;
      }
    } else {
      real_x = x+static_cast<int>(round(image_slope_x*static_cast<double>(x)))+image_offset_x;
      real_y = y+static_cast<int>(round(image_slope_y*static_cast<double>(y)))+image_offset_y;
    }

    return true;
  }

  bool inline onimage(int x,int y) const {
    int real_x=0;
    int real_y=0;

    to_real(x,y,real_x,real_y);

    if((real_x >= 0) &&
       (real_y >= 0) &&
       (real_x <  image_width()) &&
       (real_y <  image_height())) {
      return true;
    } else {
      return false;
    }
  }
  
  // Really int should be a template argument for position type
  SwiftImagePosition<int> find_image_offset(const SwiftImage<_prec> &in, // Compare to this image
                                            int threshold,               // At most this distance away
                                            int increment=1,             // Use every increment pixels (speeds things up)
                                            int useblock=0               // See multiply_sum_shift
                                           ) const {
 
    if(increment==0) {
      throw (std::invalid_argument("ERROR: Increment to find_image_offset is 0 will loop forever"));
    }

    double sum=0;
    double maxsum=0;
    int    maxsum_x_offset=0;
    int    maxsum_y_offset=0;

    bool first=true;

    vector<vector<double> > multiply_results((threshold*2)+1,vector<double>((threshold*2)+1,0));
    
#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for(int x_offset=0-threshold;x_offset<threshold;x_offset++) {
      for(int y_offset=0-threshold;y_offset<threshold;y_offset++) {
        multiply_results[x_offset+threshold][y_offset+threshold] = multiply_sum_shift(in,x_offset,y_offset,increment,useblock);
      }
    }


    for(int x_offset=0-threshold;x_offset<threshold;x_offset++) {
      for(int y_offset=0-threshold;y_offset<threshold;y_offset++) {
        
        sum = multiply_results[x_offset+threshold][y_offset+threshold];
        
        if((sum > maxsum) || first) {
          // err << "offset: " << x_offset << "," << y_offset << " sum: " << maxsum << endl;
          maxsum=sum;
          maxsum_x_offset = x_offset;
          maxsum_y_offset = y_offset;
          first=false;
        }
      }
    }

    return SwiftImagePosition<int>(maxsum_x_offset,maxsum_y_offset);
  }

  void apply_offset(SwiftImagePosition<int> offset) {
    image_offset_x += offset.x;
    image_offset_y += offset.y;
  }
  
  void apply_offset_slope(const SwiftImagePosition<int>    &offset,
                          const SwiftImagePosition<double> &slope
                         ) {
    image_offset_x += offset.x;
    image_offset_y += offset.y;

    image_slope_x += slope.x;
    image_slope_y += slope.y;
  }

  int sum() const {
  
    int this_sum=0;

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        this_sum += (*this)(x,y);
      }
    }
  
    return this_sum;
  }

  double multiply_sum_shift(const SwiftImage<_prec> &in, // Multiple against this image
                            int x_offset,                // Apply this offset (x)
                            int y_offset,                // Apply this offset (y)
                            int increment=1,             // Increment by this number of pixels across image (speed up)
                            int useblock=0               // Use a block of the image this big (in the centre)
                           ) const {

    double sum=0;


    if((useblock>image_width()) && (useblock>image_height())) useblock=0;

    _prec maxval = max();
    _prec inmax  = in.max();
    if(inmax > maxval) maxval = inmax; 

    int x_start=min_x();
    int y_start=min_y();
    int x_max = max_x();
    int y_max = max_y();

    if(useblock != 0) {
      x_start = (image_width()/2)-(useblock/2);
      y_start = (image_height()/2)-(useblock/2);
      x_max = x_start+useblock;
      y_max = y_start+useblock;
      //cout << "x_start: " << x_start << endl;
      //cout << "x_max  : " << x_max   << endl;
    }

    for(int x=x_start;x<x_max;x+=increment) {
      for(int y=y_start;y<y_max;y+=increment) {
        if(onimage(x,y) && in.onimage(x+x_offset,y+y_offset)) {
          sum += (static_cast<double>((*this)(x,y)) * static_cast<double>(in(x+x_offset,y+y_offset)))/maxval;
        }
      }
    }

    return sum;
  }

  SwiftImage<_prec> operator*(const SwiftImage<_prec> &rhs) const {
    SwiftImage<_prec> ret(0,0);

    ret=(*this);

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        if(rhs.onimage(x,y)) {
          ret(x,y) = (*this)(x,y) * rhs(x,y);
        } else ret(x,y)=0;
      }
    }

    return ret;
  }
  
  SwiftImage<_prec> operator*(int rhs) const {
    SwiftImage<_prec> ret(image_width(),image_height());
    ret = (*this);
    
    for(size_t n=0;n<ret.image.size();n++) {
     ret.image[n] = ret.image[n]*rhs;
    }
    
    return ret;
  }

  SwiftImage<_prec> dilate() {
    
    SwiftImage<_prec> img = (*this);

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        bool addone=false;

        for(int cx=-1;cx<=1;cx++) {
          for(int cy=-1;cy<=1;cy++) {
            if(img.onimage(x+cx,y+cy)) if((*this)(x+cx,y+cy) > 0) addone=true;
          }
        }

        if(addone) {
          img(x,y) = img(x,y) + 1;
        }
      }
    }

    return img;
  }


  inline _prec max() const {
   
    _prec max=image[0];

    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        _prec c = (*this)(x,y);
        if(c > max) max = c;
      }
    }

    cache_max = max;
    cache_max_ok=false;// reenable max caching at some point

    return max;
  }

  int min_x() const {
    return 0-image_offset_x;
  }

  int max_x() const {
    return image_width()-image_offset_x;
  }

  int min_y() const {
    return 0-image_offset_y;
  }

  int max_y() const {
    return image_height()-image_offset_y;
  }

  SwiftImage<_prec> subimage_offset(int x_min,
                                    int x_max,
                                    int y_min,
                                    int y_max,
                                    int offset_x,
                                    int offset_y) {
    
    SwiftImage<_prec> offseted = (*this);

    for(int x=x_min;x<x_max;x++) {
      for(int y=y_min;y<y_max;y++) {
        if(onimage(x+offset_x,y+offset_y)) {
          offseted(x,y) = (*this)(x+offset_x,y+offset_y);
        } else {
          offseted(x,y) = 0;
        }
      }
    }

    return offseted;
  }
  
  SwiftImage<_prec> crop(int x_min,int x_max,int y_min,int y_max) const {
    
    int x_size = x_max-x_min;
    int y_size = y_max-y_min;
    SwiftImage<_prec> cropimg(x_size,y_size);
    cropimg.apply_offset(SwiftImagePosition<>(0-x_min,0-y_min));
    
    for(int x=x_min;x<x_max;x++) {
      for(int y=y_min;y<y_max;y++) {
	cropimg(x,y) = (*this)(x,y);
      }
    }
    
    return cropimg;
  }

  SwiftImage<_prec> doublesize() {
    image_offset_x=0;
    image_offset_y=0;
    
    SwiftImage<_prec> dimg(image_width()*2,image_height()*2);
    
    for(int x=min_x();x<max_x();x++) {
      for(int y=min_y();y<max_y();y++) {
        dimg(x*2    ,y*2    ) = (*this)(x,y);
        dimg((x*2)+1,y*2    ) = (*this)(x,y);
        dimg(x*2    ,(y*2)+1) = (*this)(x,y);
        dimg((x*2)+1,(y*2)+1) = (*this)(x,y);
      }
    }
    
    return dimg;
  }

  void apply_offset_map(const vector<vector<SwiftImagePosition<> > > &offset_map_in) {
    if(offset_map_enable==false) { 
      offset_map_enable=true;
      offset_map = offset_map_in;
    } else {
      //TODO: If the offset maps are of different sizes this will fail completely
      if(offset_map.size() != offset_map_in.size()) cerr << "************************************ ERROR OFFSET MAPS NOT OF SAME SIZE" << endl;
      for(size_t x=0;x<offset_map.size();x++) {
        for(size_t y=0;y<offset_map[x].size();y++) {
          offset_map[x][y] += offset_map_in[x][y];
        }
      }
    }
  }

  // should really return a size_t
  inline int image_width() const {
    return m_image_width;
  }

  inline int image_height() const {
    return m_image_height;
  }

  inline void clear_offset() {
    image_offset_x=0;
    image_offset_y=0;
    offset_map_enable=false;
    offset_map.clear();
  }

  inline size_t get_index(int x,int y) const {
    int real_x = x+image_offset_x;
    int real_y = y+image_offset_y;

    return (real_y*m_image_width)+real_x;
  }

  inline size_t max_index() const {
    return image.size();
  }

  inline SwiftImagePosition<> get_offset() const {
    return SwiftImagePosition<>(image_offset_x,image_offset_y);
  }

  void dump_offset_map(ostream &out) {

    out << "X MAP:" << endl;
    for(size_t x=0;x<offset_map.size();x++) {
      for(size_t y=0;y<offset_map[x].size();y++) {
        out << setw(3) << offset_map[x][y].x << " ";
      }
      out << endl;
    }
    
    out << "Y MAP:" << endl;
    for(size_t x=0;x<offset_map.size();x++) {
      for(size_t y=0;y<offset_map[x].size();y++) {
        out << setw(3) << offset_map[x][y].y << " ";
      }
      out << endl;
    }
  }
  
  string dump_offset_map_xml() {

    ostringstream str;

    str << "<xmap>" << endl;
    for(size_t x=0;x<offset_map.size();x++) {
      for(size_t y=0;y<offset_map[x].size();y++) {
        str << "<offset x=\"" << x << "\" y=\"" << y << "\" value=\"" << offset_map[x][y].x << "\" />";
      }
      str << endl;
    }
    str << "</xmap>" << endl;
    
    str << "<ymap>" << endl;
    for(size_t x=0;x<offset_map.size();x++) {
      for(size_t y=0;y<offset_map[x].size();y++) {
        str << "<offset x=\"" << x << "\" y=\"" << y << "\" value=\"" << offset_map[x][y].x << "\" />";
      }
      str << endl;
    }
    str << "</ymap>" << endl;

    return str.str();
  }

  bool dump(ostream &out) const;

  // Fix protection, many of these can be made protected (this should be friend of itself)

  vector<_prec> image;
  unsigned int m_image_width;
  unsigned int m_image_height;
  int image_offset_x;
  int image_offset_y;

  double image_slope_x;
  double image_slope_y;

  mutable bool cache_max_ok;
  mutable _prec cache_max;
  Timetagger m_tt;
  vector<vector<SwiftImagePosition<> > > offset_map;
  bool offset_map_enable;
  ostream &err;
};

#include "SwiftImage.cpp"

#endif
