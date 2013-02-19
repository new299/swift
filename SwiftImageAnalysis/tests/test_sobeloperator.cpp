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
#include "test_sobeloperator.h"
#include "SwiftImage.h"
#include "SobelOperator.h"
#include "MorphologicalOpening.h"
#include "MorphologicalClosing.h"
#include "AdaptiveThreshold.h"
#include "NWThreshold.h"

#include <iostream>

void test_sobeloperator(UnitTest &ut) {

  ut.begin_test_set("SobelOperator");

  SwiftImage<uint16> img("./Images/tiny5dot.tif");
  SobelOperator<uint16> so;
  SwiftImage<uint16> imgso = so.process(img);
  ut.test(imgso(0,0),static_cast<uint16>(0));
  ut.test(imgso(1,0),static_cast<uint16>(0));
  ut.test(imgso(2,0),static_cast<uint16>(0));
  ut.test(imgso(3,0),static_cast<uint16>(0));
  ut.test(imgso(4,0),static_cast<uint16>(0));

  ut.test(imgso(0,1),static_cast<uint16>(0));
  ut.test(imgso(1,1),static_cast<uint16>(65534));
  ut.test(imgso(2,1),static_cast<uint16>(65534));
  ut.test(imgso(3,1),static_cast<uint16>(65534));
  ut.test(imgso(4,1),static_cast<uint16>(0));

  ut.test(imgso(0,2),static_cast<uint16>(0));
  ut.test(imgso(1,2),static_cast<uint16>(65534));
  ut.test(imgso(2,2),static_cast<uint16>(0));
  ut.test(imgso(3,2),static_cast<uint16>(65534));
  ut.test(imgso(4,2),static_cast<uint16>(0));

  ut.test(imgso(0,3),static_cast<uint16>(0));
  ut.test(imgso(1,3),static_cast<uint16>(65534));
  ut.test(imgso(2,3),static_cast<uint16>(65534));
  ut.test(imgso(3,3),static_cast<uint16>(65534));
  ut.test(imgso(4,3),static_cast<uint16>(0));

  ut.test(imgso(0,4),static_cast<uint16>(0));
  ut.test(imgso(1,4),static_cast<uint16>(0));
  ut.test(imgso(2,4),static_cast<uint16>(0));
  ut.test(imgso(3,4),static_cast<uint16>(0));
  ut.test(imgso(4,4),static_cast<uint16>(0));
/*
  SwiftImage<uint16> i2("./Images/run475_lane1tile1/C2.1/s_1_1_a.tif");
  MorphologicalClosing<uint16> mc(4);
  SwiftImage<uint16> i4 = mc.process(i2);
 
  i2 = i2/((i4/10000)+1);
  SwiftImage<uint16> i3 = so.process(i2);
  
  i3.save("./sobel_morph.tif");

  NWThreshold<uint16> at(2,10);
  i3 = at.process(i3);
  i3.save("./sobel_morph_adaptive.tif");
*/
  ut.end_test_set();
}

