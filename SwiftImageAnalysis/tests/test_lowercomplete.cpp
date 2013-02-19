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
#include "test_lowercomplete.h"
#include "SwiftImage.h"
#include "LowerComplete.h"
#include "NWThreshold.h"

#include <iostream>

void test_lowercomplete(UnitTest &ut) {

  ut.begin_test_set("LowerCompletion");
  
  LowerComplete<uint16> lowercomp;
  
  SwiftImage<uint16> i50("./Images/5with0.tif");
  SwiftImage<uint16> i50_lowercomp = lowercomp.process(i50);
  
  SwiftImage<uint16> i11("./Images/tiny11.tif");
  
  SwiftImage<uint16> i11_lowercomp = lowercomp.process(i11);

  SwiftImage<uint16> i13("./Images/tiny13.tif");
  SwiftImage<uint16> i13_lowercomp = lowercomp.process(i13);

  SwiftImage<uint16> il("./Images/tinyline.tif");
  SwiftImage<uint16> il_lowercomp = lowercomp.process(il);

  ut.end_test_set();
}

