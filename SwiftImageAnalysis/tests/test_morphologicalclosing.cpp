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
#include "test_morphologicalclosing.h"
#include "SwiftImage.h"
#include "MorphologicalClosing.h"

#include <iostream>

void test_morphologicalclosing(UnitTest &ut) {

  ut.begin_test_set("MorphologicalClosing");
/*
  MorphologicalClosing<uint16> mo3(2,false,MorphologicalClosing<uint16>::mask_type_circle);
  SwiftImage<uint16> i2("./Images/run475_lane1tile1/C2.1/s_1_1_a.tif");
  SwiftImage<uint16> i3 = mo3.process(i2);

  SwiftImage<uint16> i4 = i2/((i3/10000)+1);
  i2.save("./premorphc.tif");
  i3.save("./morphc.tif");
  i4.save("./divmorphc.tif");
*/
  

  ut.end_test_set();
}

