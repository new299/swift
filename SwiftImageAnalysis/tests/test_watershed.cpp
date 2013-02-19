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
#include "test_watershed.h"
#include "SwiftImage.h"
#include "EuclideanDistanceMap.h"
#include "Watershed.h"
#include "NWThreshold.h"
#include "Invert.h"

#include <iostream>

void test_watershed(UnitTest &ut) {

  ut.begin_test_set("Watershed");
  
  EuclideanDistanceMap<uint16> edm;
  Watershed<uint16> wat;
  
  SwiftImage<uint16> small_frompaper("./Images/watsmall.tif");
  SwiftImage<uint16> sfp_wat = wat.process(small_frompaper);
  ut.test(sfp_wat(0,0),static_cast<uint16>(1));
  ut.test(sfp_wat(1,0),static_cast<uint16>(1));
  ut.test(sfp_wat(2,0),static_cast<uint16>(0));
  ut.test(sfp_wat(3,0),static_cast<uint16>(1));
  ut.test(sfp_wat(4,0),static_cast<uint16>(1));
  
  ut.test(sfp_wat(0,1),static_cast<uint16>(1));
  ut.test(sfp_wat(1,1),static_cast<uint16>(1));
  ut.test(sfp_wat(2,1),static_cast<uint16>(0));
  ut.test(sfp_wat(3,1),static_cast<uint16>(1));
  ut.test(sfp_wat(4,1),static_cast<uint16>(1));
  
  ut.test(sfp_wat(0,2),static_cast<uint16>(0));
  ut.test(sfp_wat(1,2),static_cast<uint16>(0));
  ut.test(sfp_wat(2,2),static_cast<uint16>(0));
  ut.test(sfp_wat(3,2),static_cast<uint16>(0));
  ut.test(sfp_wat(4,2),static_cast<uint16>(0));
  
  ut.test(sfp_wat(0,3),static_cast<uint16>(1));
  ut.test(sfp_wat(1,3),static_cast<uint16>(1));
  ut.test(sfp_wat(2,3),static_cast<uint16>(0));
  ut.test(sfp_wat(3,3),static_cast<uint16>(1));
  ut.test(sfp_wat(4,3),static_cast<uint16>(1));
  
  ut.test(sfp_wat(0,4),static_cast<uint16>(1));
  ut.test(sfp_wat(1,4),static_cast<uint16>(1));
  ut.test(sfp_wat(2,4),static_cast<uint16>(0));
  ut.test(sfp_wat(3,4),static_cast<uint16>(1));
  ut.test(sfp_wat(4,4),static_cast<uint16>(1));
  
  SwiftImage<uint16> i11("./Images/tiny11.tif");
  SwiftImage<uint16> i11_edm = edm.process(i11);
  
  
  Invert<uint16> inv;
  SwiftImage<uint16> i11_edm_inv = inv.process(i11_edm);
  
  SwiftImage<uint16> i11_wat = wat.process(i11_edm_inv);
 
  SwiftImage<uint16> f = i11 && i11_wat;

  ut.test(i11_wat(0,0 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,1 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,9 ),static_cast<uint16>(1));
  ut.test(i11_wat(0,10),static_cast<uint16>(1));
  
  ut.test(i11_wat(1,0 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,1 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,9 ),static_cast<uint16>(1));
  ut.test(i11_wat(1,10),static_cast<uint16>(1));
  
  ut.test(i11_wat(2,0 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,1 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,9 ),static_cast<uint16>(1));
  ut.test(i11_wat(2,10),static_cast<uint16>(1));
  
  ut.test(i11_wat(3,0 ),static_cast<uint16>(0));
  ut.test(i11_wat(3,1 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,9 ),static_cast<uint16>(1));
  ut.test(i11_wat(3,10),static_cast<uint16>(0));
  
  ut.test(i11_wat(4,0 ),static_cast<uint16>(0));
  ut.test(i11_wat(4,1 ),static_cast<uint16>(0));
  ut.test(i11_wat(4,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(4,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(4,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(4,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(4,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(4,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(4,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(4,9 ),static_cast<uint16>(0));
  ut.test(i11_wat(4,10),static_cast<uint16>(0));
  
  ut.test(i11_wat(5,0 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,1 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,2 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,3 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,4 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,5 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,6 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,7 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,8 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,9 ),static_cast<uint16>(0));
  ut.test(i11_wat(5,10),static_cast<uint16>(0));
  
  ut.test(i11_wat(6,0 ),static_cast<uint16>(0));
  ut.test(i11_wat(6,1 ),static_cast<uint16>(0));
  ut.test(i11_wat(6,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(6,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(6,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(6,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(6,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(6,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(6,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(6,9 ),static_cast<uint16>(0));
  ut.test(i11_wat(6,10),static_cast<uint16>(0));
  
  ut.test(i11_wat(7,0 ),static_cast<uint16>(0));
  ut.test(i11_wat(7,1 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,9 ),static_cast<uint16>(1));
  ut.test(i11_wat(7,10),static_cast<uint16>(0));
  
  ut.test(i11_wat(8,0 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,1 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,9 ),static_cast<uint16>(1));
  ut.test(i11_wat(8,10),static_cast<uint16>(1));
  
  ut.test(i11_wat(9,0 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,1 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,9 ),static_cast<uint16>(1));
  ut.test(i11_wat(9,10),static_cast<uint16>(1));
  
  ut.test(i11_wat(10,0 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,1 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,2 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,3 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,4 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,5 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,6 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,7 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,8 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,9 ),static_cast<uint16>(1));
  ut.test(i11_wat(10,10),static_cast<uint16>(1));
  
  SwiftImage<uint16> i13("./Images/tiny13.tif");

  SwiftImage<uint16> i13_edm = edm.process(i13);
  SwiftImage<uint16> i13_inv = inv.process(i13_edm);
  SwiftImage<uint16> i13_wat = wat.process(i13_inv);
  SwiftImage<uint16> i13_cmb = i13 && i13_wat;
  i13_cmb.save("wat13_cmb.tif");
 
/*
  NWThreshold<uint16> nwt(5,0.5,2000,NWThreshold<uint16>::mask_type_square);//0.75 beautiful

  SwiftImage<uint16> i1("./Images/run475_lane1tile1/C2.1/s_1_1_a.tif");
  SwiftImage<uint16> i2 = nwt.process(i1);
  SwiftImage<uint16> i3 = edm.process(i2);
  SwiftImage<uint16> i4 = inv.process(i3);
  SwiftImage<uint16> i5 = wat.process(i4);
  SwiftImage<uint16> i6 = i2 && i5;
  
  i2.save("nwt_watershed.tif");
  i4.save("inv_watershed.tif");
  i5.save("wat_watershed.tif");
  i6.save("watershed.tif");
*/
  ut.end_test_set();
}
