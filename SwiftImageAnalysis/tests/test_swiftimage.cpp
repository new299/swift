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

#include "utf.h"
#include "test_swiftimage.h"
#include "SwiftImage.h"

#include <iostream>

void test_swiftimage(UnitTest &ut) {

  ut.begin_test_set("SwiftImage");
  SwiftImage<uint16> i("./Images/tiny2.tif");
  // i.dump(std::cout);

  ut.test(i(0,0),static_cast<uint16>(0));
  ut.test(i(1,0),static_cast<uint16>(0));
  ut.test(i(0,1),static_cast<uint16>(65535));
  ut.test(i(1,1),static_cast<uint16>(65535));

  SwiftImage<uint16> id2("./Images/tinyd2.tif");

  i = i/id2;
  ut.test(i(0,0),static_cast<uint16>(0));
  ut.test(i(1,0),static_cast<uint16>(0));
  ut.test(i(0,1),static_cast<uint16>(32767));
  ut.test(i(1,1),static_cast<uint16>(32767));
  
  SwiftImage<uint16> i2("./Images/tiny4.tif");
  ut.test(i2(0,0),static_cast<uint16>(0));
  ut.test(i2(1,0),static_cast<uint16>(0));
  ut.test(i2(2,0),static_cast<uint16>(0));
  ut.test(i2(3,0),static_cast<uint16>(0));
  ut.test(i2(0,1),static_cast<uint16>(0));
  ut.test(i2(1,1),static_cast<uint16>(0));
  ut.test(i2(2,1),static_cast<uint16>(0));
  ut.test(i2(3,1),static_cast<uint16>(0));
  ut.test(i2(0,2),static_cast<uint16>(65535));
  ut.test(i2(1,2),static_cast<uint16>(65535));
  ut.test(i2(2,2),static_cast<uint16>(65535));
  ut.test(i2(3,2),static_cast<uint16>(65535));
  ut.test(i2(0,3),static_cast<uint16>(65535));
  ut.test(i2(1,3),static_cast<uint16>(65535));
  ut.test(i2(2,3),static_cast<uint16>(65535));
  ut.test(i2(3,3),static_cast<uint16>(65535));


  ut.end_test_set();

  test_swiftimage_find_image_offset(ut);
}

void test_swiftimage_find_image_offset(UnitTest &ut) {

  ut.begin_test_set("SwiftImage - find_image_offset");
  
  // Test offseting
  SwiftImage<uint16> i3("./Images/dot5.tif");
  SwiftImage<uint16> i4("./Images/dot5shift1.tif");

  SwiftImagePosition<> offset = i3.find_image_offset(i4,4,1,0);
  ut.test(offset.x,1);
  ut.test(offset.y,-1);
  
  //SwiftImagePosition<> offset_fft = i3.find_image_offset_fft(i4);
  //ut.test(offset_fft.x,1);
  //ut.test(offset_fft.y,-1);

  ut.end_test_set();
}
