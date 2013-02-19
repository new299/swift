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
#include "test_segmentation.h"
#include "SwiftImage.h"
#include "Segmentation.h"
#include "RunLengthEncode.h"
#include "EuclideanDistanceMap.h"
#include "RLERun.h"
#include "Watershed.h"
#include "NWThreshold.h"
#include "Invert.h"

#include <iostream>

void test_segmentation(UnitTest &ut) {

  ut.begin_test_set("Segmentation");

  SwiftImage<uint16> img("./Images/tinyline.tif");
  
  Segmentation<uint16> segmenter;

  vector<SwiftImageObject<> > objs = segmenter.process(img);

  ut.test(static_cast<int>(objs.size()),1);
  ut.test(objs[0].pixels[0].pos.x,2);
  ut.test(objs[0].pixels[0].pos.y,1);
  ut.test(objs[0].pixels[0].length,9);

  /*
  EuclideanDistanceMap<uint16> edm;
  Watershed<uint16> wat;
  Invert<uint16> inv;
  NWThreshold<uint16> nwt(5,0.5,60000,NWThreshold<uint16>::mask_type_square);//0.75 beautiful

  SwiftImage<uint16> i1("./Images/run475_lane1tile1/C2.1/s_1_1_a.tif");
  cout << "Load complete" << endl;
  SwiftImage<uint16> i2 = nwt.process(i1);
  cout << "Thresholding complete" << endl;
  i2.save("./nwt.tiff");
  SwiftImage<uint16> i3 = edm.process(i2);
  cout << "EDM complete" << endl;
  SwiftImage<uint16> i4 = inv.process(i3);
  cout << "INV complete" << endl;
  SwiftImage<uint16> i5 = wat.process(i4);
  cout << "Wat complete" << endl;
  SwiftImage<uint16> i6 = i2 && i5;

  vector<SwiftImageObject<> > objs1 = segmenter.process(i6);
  ut.test(static_cast<int>(objs1.size()),1);
  cout << "Seg complete" << endl;
  */


  ut.end_test_set();
}

