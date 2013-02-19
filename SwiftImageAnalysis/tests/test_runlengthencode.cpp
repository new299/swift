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
#include "test_runlengthencode.h"
#include "SwiftImage.h"
#include "RunLengthEncode.h"
#include "RLERun.h"

#include <iostream>

void test_runlengthencode(UnitTest &ut) {

  ut.begin_test_set("RunLengthEncode");

  SwiftImage<uint16> img("./Images/tinyline.tif");
  RunLengthEncode<uint16> rle;

  vector<RLERun<> > runs = rle.process(img);
 
  ut.test(static_cast<int>(runs.size())   ,1);
  ut.test(runs[0].pos.x ,2);
  ut.test(runs[0].pos.y ,1);
  ut.test(runs[0].length,9);

  ut.end_test_set();
}

