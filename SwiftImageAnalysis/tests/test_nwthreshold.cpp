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
#include "test_nwthreshold.h"
#include "SwiftImage.h"
#include "NWThreshold.h"
#include "AdaptiveThreshold.h"
#include "MorphologicalOpening.h"
#include "MorphologicalClosing.h"

#include <iostream>

void test_nwthreshold(UnitTest &ut) {

  ut.begin_test_set("NWThreshold");

  SwiftImage<uint16> img("./Images/tiny5dot.tif");
  NWThreshold<uint16> nwt(5,0.5,65534,NWThreshold<uint16>::mask_type_square);//0.75 beautiful
  SwiftImage<uint16> imgat = nwt.process(img);
  imgat.save("./nwthreshold_tiny5dot.tif");

  // AdaptiveThreshold<uint16> at(5,1);

/*  SwiftImage<uint16> i1("./Images/run475_lane1tile1/C2.1/s_1_1_a.tif");
  // MorphologicalClosing<uint16> mc(4);
  // SwiftImage<uint16> i4 = mc.process(i1);
 
  // SwiftImage<uint16> i2 = i1/((i4/10000)+1);
  // SwiftImage<uint16> i3 = at.process(i2);
  
  SwiftImage<uint16> i5 = nwt.process(i1);

  // i3.save("./nwthreshold_morph.tif");
  i5.save("./nwthreshold.tif");
*/
  ut.end_test_set();
}

