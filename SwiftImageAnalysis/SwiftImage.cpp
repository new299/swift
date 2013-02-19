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

using namespace std;

/// This method saves a tiff image
template<class _prec>
bool SwiftImage<_prec>::save(const char *filename) const {

  SwiftImage<_prec> offseted = (*this);
  offseted.copy_with_offset((*this));
  
  TIFF* tif = TIFFOpen(filename, "w");

  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, m_image_width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, m_image_height);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, m_image_height);

  // TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_CCITTFAX4);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

  TIFFSetField(tif, TIFFTAG_XRESOLUTION, 150.0);
  TIFFSetField(tif, TIFFTAG_YRESOLUTION, 150.0);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
  
  // Write the information to the file
  for(unsigned int n=0;n<offseted.m_image_height;n++) {
    TIFFWriteScanline(tif, &(offseted.image[offseted.m_image_width*n]),n, offseted.image.size());
  }

  TIFFClose(tif);

  return true;
}


/// This method loads a tiff image in to the image vector
template<class _prec>
bool SwiftImage<_prec>::load(const char *filename) {
  
  TIFF *tif = TIFFOpen(filename,"r");
  if(tif == NULL) return false;

  // Determine image width and height
  unsigned int width;
  unsigned int height;

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
  if(width == 0) return false;

  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
  if(height == 0) return false;

  // Clear image
  image.clear();

  long result;

  uint32 imagelength;
  uint16 *buf;

  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);

  unsigned int scanlinesize = TIFFScanlineSize(tif);
  buf = static_cast<uint16 *>(_TIFFmalloc(TIFFScanlineSize(tif)));
  
  // Read Image data
  for (uint32 row = 0; row < imagelength; row++) {
    result = TIFFReadScanline(tif,buf,row);

    uint16 *bp = buf;
    for(unsigned int n=0;n<(scanlinesize/sizeof(uint16));n++) {
      uint16 val = *(static_cast<uint16 *>(bp));
      image.push_back(val);
      bp++;
    }
    if(result == -1) {return false;}
  }
  
  m_image_width  = width;// imagelength;
  m_image_height = height;//scanlinesize/sizeof(uint16);

  _TIFFfree(buf);
    
  TIFFClose(tif);
  
  return true;
}

/// Dump image to ostream
template<class _prec>
bool SwiftImage<_prec>::dump(ostream &out) const {
  for(int y=-20;y<20;y++) {
    for(int x=-20;x<20;x++) {
      if(!onimage(x,y)) out << " 0";
                   else out << " " << (*this)(x,y);
    }
    out << endl;
  }
  
  return true;
}
