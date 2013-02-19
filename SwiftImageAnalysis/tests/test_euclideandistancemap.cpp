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
#include "test_euclideandistancemap.h"
#include "SwiftImage.h"
#include "EuclideanDistanceMap.h"
#include "NWThreshold.h"

#include <iostream>

void test_euclideandistancemap(UnitTest &ut) {

  ut.begin_test_set("EuclideanDistanceMap");

  EuclideanDistanceMap<uint16> edm;
  
  SwiftImage<uint16> i10("./Images/tiny10.tif");
  SwiftImage<uint16> i10_edm = edm.process(i10);
  
  ut.test(i10_edm(0,0),static_cast<uint16>(0));
  ut.test(i10_edm(0,1),static_cast<uint16>(0));
  ut.test(i10_edm(0,2),static_cast<uint16>(0));
  ut.test(i10_edm(0,3),static_cast<uint16>(0));
  ut.test(i10_edm(0,4),static_cast<uint16>(0));
  ut.test(i10_edm(0,5),static_cast<uint16>(0));
  ut.test(i10_edm(0,6),static_cast<uint16>(0));
  ut.test(i10_edm(0,7),static_cast<uint16>(0));
  ut.test(i10_edm(0,8),static_cast<uint16>(0));
  ut.test(i10_edm(0,9),static_cast<uint16>(0));
  
  ut.test(i10_edm(1,0),static_cast<uint16>(0));
  ut.test(i10_edm(1,1),static_cast<uint16>(1));
  ut.test(i10_edm(1,2),static_cast<uint16>(1));
  ut.test(i10_edm(1,3),static_cast<uint16>(1));
  ut.test(i10_edm(1,4),static_cast<uint16>(1));
  ut.test(i10_edm(1,5),static_cast<uint16>(1));
  ut.test(i10_edm(1,6),static_cast<uint16>(1));
  ut.test(i10_edm(1,7),static_cast<uint16>(1));
  ut.test(i10_edm(1,8),static_cast<uint16>(1));
  ut.test(i10_edm(1,9),static_cast<uint16>(0));
  
  ut.test(i10_edm(2,0),static_cast<uint16>(0));
  ut.test(i10_edm(2,1),static_cast<uint16>(1));
  ut.test(i10_edm(2,2),static_cast<uint16>(2));
  ut.test(i10_edm(2,3),static_cast<uint16>(2));
  ut.test(i10_edm(2,4),static_cast<uint16>(2));
  ut.test(i10_edm(2,5),static_cast<uint16>(2));
  ut.test(i10_edm(2,6),static_cast<uint16>(2));
  ut.test(i10_edm(2,7),static_cast<uint16>(2));
  ut.test(i10_edm(2,8),static_cast<uint16>(1));
  ut.test(i10_edm(2,9),static_cast<uint16>(0));
  
  ut.test(i10_edm(3,0),static_cast<uint16>(0));
  ut.test(i10_edm(3,1),static_cast<uint16>(1));
  ut.test(i10_edm(3,2),static_cast<uint16>(2));
  ut.test(i10_edm(3,3),static_cast<uint16>(3));
  ut.test(i10_edm(3,4),static_cast<uint16>(3));
  ut.test(i10_edm(3,5),static_cast<uint16>(3));
  ut.test(i10_edm(3,6),static_cast<uint16>(3));
  ut.test(i10_edm(3,7),static_cast<uint16>(2));
  ut.test(i10_edm(3,8),static_cast<uint16>(1));
  ut.test(i10_edm(3,9),static_cast<uint16>(0));
  
  ut.test(i10_edm(4,0),static_cast<uint16>(0));
  ut.test(i10_edm(4,1),static_cast<uint16>(1));
  ut.test(i10_edm(4,2),static_cast<uint16>(2));
  ut.test(i10_edm(4,3),static_cast<uint16>(3));
  ut.test(i10_edm(4,4),static_cast<uint16>(4));
  ut.test(i10_edm(4,5),static_cast<uint16>(4));
  ut.test(i10_edm(4,6),static_cast<uint16>(3));
  ut.test(i10_edm(4,7),static_cast<uint16>(2));
  ut.test(i10_edm(4,8),static_cast<uint16>(1));
  ut.test(i10_edm(4,9),static_cast<uint16>(0));
 
  ut.test(i10_edm(5,0),static_cast<uint16>(0));
  ut.test(i10_edm(5,1),static_cast<uint16>(1));
  ut.test(i10_edm(5,2),static_cast<uint16>(2));
  ut.test(i10_edm(5,3),static_cast<uint16>(3));
  ut.test(i10_edm(5,4),static_cast<uint16>(4));
  ut.test(i10_edm(5,5),static_cast<uint16>(4));
  ut.test(i10_edm(5,6),static_cast<uint16>(3));
  ut.test(i10_edm(5,7),static_cast<uint16>(2));
  ut.test(i10_edm(5,8),static_cast<uint16>(1));
  ut.test(i10_edm(5,9),static_cast<uint16>(0));
  
  ut.test(i10_edm(6,0),static_cast<uint16>(0));
  ut.test(i10_edm(6,1),static_cast<uint16>(1));
  ut.test(i10_edm(6,2),static_cast<uint16>(2));
  ut.test(i10_edm(6,3),static_cast<uint16>(3));
  ut.test(i10_edm(6,4),static_cast<uint16>(3));
  ut.test(i10_edm(6,5),static_cast<uint16>(3));
  ut.test(i10_edm(6,6),static_cast<uint16>(3));
  ut.test(i10_edm(6,7),static_cast<uint16>(2));
  ut.test(i10_edm(6,8),static_cast<uint16>(1));
  ut.test(i10_edm(6,9),static_cast<uint16>(0));
  
  ut.test(i10_edm(7,0),static_cast<uint16>(0));
  ut.test(i10_edm(7,1),static_cast<uint16>(1));
  ut.test(i10_edm(7,2),static_cast<uint16>(2));
  ut.test(i10_edm(7,3),static_cast<uint16>(2));
  ut.test(i10_edm(7,4),static_cast<uint16>(2));
  ut.test(i10_edm(7,5),static_cast<uint16>(2));
  ut.test(i10_edm(7,6),static_cast<uint16>(2));
  ut.test(i10_edm(7,7),static_cast<uint16>(2));
  ut.test(i10_edm(7,8),static_cast<uint16>(1));
  ut.test(i10_edm(7,9),static_cast<uint16>(0));
  
  ut.test(i10_edm(8,0),static_cast<uint16>(0));
  ut.test(i10_edm(8,1),static_cast<uint16>(1));
  ut.test(i10_edm(8,2),static_cast<uint16>(1));
  ut.test(i10_edm(8,3),static_cast<uint16>(1));
  ut.test(i10_edm(8,4),static_cast<uint16>(1));
  ut.test(i10_edm(8,5),static_cast<uint16>(1));
  ut.test(i10_edm(8,6),static_cast<uint16>(1));
  ut.test(i10_edm(8,7),static_cast<uint16>(1));
  ut.test(i10_edm(8,8),static_cast<uint16>(1));
  ut.test(i10_edm(8,9),static_cast<uint16>(0));
  
  ut.test(i10_edm(9,0),static_cast<uint16>(0));
  ut.test(i10_edm(9,1),static_cast<uint16>(0));
  ut.test(i10_edm(9,2),static_cast<uint16>(0));
  ut.test(i10_edm(9,3),static_cast<uint16>(0));
  ut.test(i10_edm(9,4),static_cast<uint16>(0));
  ut.test(i10_edm(9,5),static_cast<uint16>(0));
  ut.test(i10_edm(9,6),static_cast<uint16>(0));
  ut.test(i10_edm(9,7),static_cast<uint16>(0));
  ut.test(i10_edm(9,8),static_cast<uint16>(0));
  ut.test(i10_edm(9,9),static_cast<uint16>(0));
  
  // Trying things out...
/*  NWThreshold<uint16> nwt(5,0.5,2000,NWThreshold<uint16>::mask_type_square);//0.75 beautiful

  SwiftImage<uint16> i1("./Images/run475_lane1tile1/C2.1/s_1_1_a.tif");
  SwiftImage<uint16> i2 = nwt.process(i1);

  SwiftImage<uint16> i3 = edm.process(i2);
  i3.save("edm.tiff");
*/
  ut.end_test_set();
}

