/*
    Swift (c) 2008 Genome Research Ltd.
    Authors: Nava Whiteford and Tom Skelly (new@sgenomics.org ts6@sanger.ac.uk)

    This file is part of Swift (http://swiftng.sourceforge.net).

    Swift is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Swift.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "test_channeloffsets.h"
#include <iostream>

void test_channeloffsets(UnitTest &ut) {

  int size=100;
  int num_cycles=10;
  ut.begin_test_set("ChannelOffsets");
  ChannelOffsets<uint16> co_fft   (ChannelOffsets<uint16>::fft);
  ChannelOffsets<uint16> co_simple(ChannelOffsets<uint16>::simple,size);

  vector<SwiftImage<uint16> > images;

  int size_x=size;
  int size_y=size;

  srand(12345);

  for(int n=0;n<num_cycles;n++) {
    SwiftImage<uint16> rnd(size_x,size_y);
    
    for(int r=0; r<2500; r++) {
      int x = rand() % size_x;
      int y = rand() % size_y;

      rnd(x,y) = rand() % 32768;
    }

    images.push_back(rnd);
  }
  
  co_fft   .process(images[0],images,num_cycles);
  co_simple.process(images[0],images,num_cycles);
  
  for(int n=0;n<co_fft.get_offsets().size();n++) {
    ut.test(co_fft.get_offsets()[n].x,co_simple.get_offsets()[n].x);
    ut.test(co_fft.get_offsets()[n].y,co_simple.get_offsets()[n].y);
  }

  ut.end_test_set();
}

