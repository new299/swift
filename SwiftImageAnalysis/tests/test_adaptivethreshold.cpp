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
#include "test_adaptivethreshold.h"
#include "SwiftImage.h"
#include "AdaptiveThreshold.h"
#include "MorphologicalOpening.h"

#include <iostream>

void test_adaptivethreshold(UnitTest &ut) {

  ut.begin_test_set("AdaptiveThreshold");
/*
  SwiftImage<uint16> img("./Images/tiny5dot.tif");
  AdaptiveThreshold<uint16> at(1,0.7,AdaptiveThreshold<uint16>::mask_type_circle);
  SwiftImage<uint16> imgat = at.process(img);
  ut.test(imgat(0,0),static_cast<uint16>(0));
  ut.test(imgat(1,0),static_cast<uint16>(0));
  ut.test(imgat(2,0),static_cast<uint16>(0));
  ut.test(imgat(3,0),static_cast<uint16>(0));
  ut.test(imgat(4,0),static_cast<uint16>(0));

  ut.test(imgat(0,1),static_cast<uint16>(0));
  ut.test(imgat(1,1),static_cast<uint16>(65534));
  ut.test(imgat(2,1),static_cast<uint16>(65534));
  ut.test(imgat(3,1),static_cast<uint16>(65534));
  ut.test(imgat(4,1),static_cast<uint16>(0));

  ut.test(imgat(0,2),static_cast<uint16>(0));
  ut.test(imgat(1,2),static_cast<uint16>(65534));
  ut.test(imgat(2,2),static_cast<uint16>(0));
  ut.test(imgat(3,2),static_cast<uint16>(65534));
  ut.test(imgat(4,2),static_cast<uint16>(0));

  ut.test(imgat(0,3),static_cast<uint16>(0));
  ut.test(imgat(1,3),static_cast<uint16>(65534));
  ut.test(imgat(2,3),static_cast<uint16>(65534));
  ut.test(imgat(3,3),static_cast<uint16>(65534));
  ut.test(imgat(4,3),static_cast<uint16>(0));

  ut.test(imgat(0,4),static_cast<uint16>(0));
  ut.test(imgat(1,4),static_cast<uint16>(0));
  ut.test(imgat(2,4),static_cast<uint16>(0));
  ut.test(imgat(3,4),static_cast<uint16>(0));
  ut.test(imgat(4,4),static_cast<uint16>(0));

  SwiftImage<uint16> i2("./Images/run475_lane1tile1/C2.1/s_1_1_a.tif");
  MorphologicalOpening<uint16> mo(4);
  SwiftImage<uint16> i4 = mo.process(i2);
 
  i2 = i2 - i4;
  SwiftImage<uint16> i3 = at.process(i2);
  

  i3.save("./adaptive_morph.tif");
*/
  ut.end_test_set();
}

