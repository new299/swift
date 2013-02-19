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
#include "test_morphologicalopening.h"
#include "SwiftImage.h"
#include "MorphologicalOpening.h"

#include <iostream>

void test_morphologicalopening(UnitTest &ut) {

  ut.begin_test_set("MorphologicalOpening");
  SwiftImage<uint16> img("./Images/tiny2.tif");
  
  MorphologicalOpening<uint16> mo(0);

  SwiftImage<uint16> iout = mo.process(img);
  ut.test(img(0,0),static_cast<uint16>(0));
  ut.test(img(1,0),static_cast<uint16>(0));
  ut.test(img(0,1),static_cast<uint16>(65535));
  ut.test(img(1,1),static_cast<uint16>(65535));

  SwiftImage<uint16> img2("./Images/tiny4.tif");
  MorphologicalOpening<uint16> mo1(1);
  SwiftImage<uint16> img2morph = mo1.process(img2);
  ut.test(img2morph(0,0),static_cast<uint16>(0));
  ut.test(img2morph(1,0),static_cast<uint16>(0));
  ut.test(img2morph(2,0),static_cast<uint16>(0));
  ut.test(img2morph(3,0),static_cast<uint16>(0));
  ut.test(img2morph(0,1),static_cast<uint16>(0));
  ut.test(img2morph(1,1),static_cast<uint16>(0));
  ut.test(img2morph(2,1),static_cast<uint16>(0));
  ut.test(img2morph(3,1),static_cast<uint16>(0));
  ut.test(img2morph(0,2),static_cast<uint16>(0));
  ut.test(img2morph(1,2),static_cast<uint16>(0));
  ut.test(img2morph(2,2),static_cast<uint16>(0));
  ut.test(img2morph(3,2),static_cast<uint16>(0));
  ut.test(img2morph(0,3),static_cast<uint16>(65535));
  ut.test(img2morph(1,3),static_cast<uint16>(65535));
  ut.test(img2morph(2,3),static_cast<uint16>(65535));
  ut.test(img2morph(3,3),static_cast<uint16>(65535));

  SwiftImage<uint16> img3("./Images/tiny5dot.tif");
  MorphologicalOpening<uint16> mo2(1);
  SwiftImage<uint16> img3morph = mo2.process(img3);
  img3morph.save("./morph5dot.tif");
  ut.test(img3morph(0,0),static_cast<uint16>(65535));
  ut.test(img3morph(1,0),static_cast<uint16>(65535));
  ut.test(img3morph(2,0),static_cast<uint16>(65535));
  ut.test(img3morph(3,0),static_cast<uint16>(65535));
  ut.test(img3morph(4,0),static_cast<uint16>(65535));

  ut.test(img3morph(0,1),static_cast<uint16>(65535));
  ut.test(img3morph(1,1),static_cast<uint16>(0));
  ut.test(img3morph(2,1),static_cast<uint16>(0));
  ut.test(img3morph(3,1),static_cast<uint16>(0));
  ut.test(img3morph(4,1),static_cast<uint16>(65535));

  ut.test(img3morph(0,2),static_cast<uint16>(65535));
  ut.test(img3morph(1,2),static_cast<uint16>(0));
  ut.test(img3morph(2,2),static_cast<uint16>(0));
  ut.test(img3morph(3,2),static_cast<uint16>(0));
  ut.test(img3morph(4,2),static_cast<uint16>(65535));

  ut.test(img3morph(0,3),static_cast<uint16>(65535));
  ut.test(img3morph(1,3),static_cast<uint16>(0));
  ut.test(img3morph(2,3),static_cast<uint16>(0));
  ut.test(img3morph(3,3),static_cast<uint16>(0));
  ut.test(img3morph(4,3),static_cast<uint16>(65535));

  ut.test(img3morph(0,4),static_cast<uint16>(65535));
  ut.test(img3morph(1,4),static_cast<uint16>(65535));
  ut.test(img3morph(2,4),static_cast<uint16>(65535));
  ut.test(img3morph(3,4),static_cast<uint16>(65535));
  ut.test(img3morph(4,4),static_cast<uint16>(65535));


/*  MorphologicalOpening<uint16> mo3(5,false,MorphologicalOpening<uint16>::mask_type_circle);
  SwiftImage<uint16> i2("./Images/run475_lane1tile1/C2.1/s_1_1_a.tif");
  SwiftImage<uint16> i3 = mo3.process(i2);

  SwiftImage<uint16> i4 = i2/((i3/10000)+1);
  i2.save("./premorph.tif");
  i3.save("./morph.tif");
  i4.save("./divmorph.tif");
*/

  ut.end_test_set();
}

