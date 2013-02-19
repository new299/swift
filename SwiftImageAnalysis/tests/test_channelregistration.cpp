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

#include <iostream>

#include "utf.h"
#include "test_channelregistration.h"
#include "SwiftImage.h"
#include "SwiftImageCluster.h"
#include "SwiftImageObject.h"
#include "Segmentation.h"
#include "ChannelRegistration.h"
#include "RunLengthEncode.h"
#include "EuclideanDistanceMap.h"
#include "RLERun.h"
#include "Watershed.h"
#include "NWThreshold.h"
#include "Invert.h"


void test_channelregistration(UnitTest &ut) {

  ut.begin_test_set("ChannelRegistration");

  Segmentation<uint16> segmenter;

  EuclideanDistanceMap<uint16> edm;
  Watershed<uint16> wat;
  Invert<uint16> inv;
  NWThreshold<uint16> nwt(5,0.5,60000,NWThreshold<uint16>::mask_type_square);//0.75 beautiful

  SwiftImage<uint16> c1_i1("./Images/smallfake2/c1_a.tif");
  cout << "Load complete" << endl;
  SwiftImage<uint16> c1_i2 = nwt.process(c1_i1);
  cout << "Thresholding complete" << endl;
  SwiftImage<uint16> c1_i3 = edm.process(c1_i2);
  cout << "EDM complete" << endl;
  SwiftImage<uint16> c1_i4 = inv.process(c1_i3);
  cout << "INV complete" << endl;
  SwiftImage<uint16> c1_i5 = wat.process(c1_i4);
  cout << "Wat complete" << endl;
  SwiftImage<uint16> c1_i6 = c1_i2 && c1_i5;
  vector<SwiftImageObject<> > c1_objs1 = segmenter.process(c1_i6);
  
  SwiftImage<uint16> c2_i1("./Images/smallfake2/c2_a.tif");
  cout << "Load complete" << endl;
  SwiftImage<uint16> c2_i2 = nwt.process(c2_i1);
  cout << "Thresholding complete" << endl;
  SwiftImage<uint16> c2_i3 = edm.process(c2_i2);
  cout << "EDM complete" << endl;
  SwiftImage<uint16> c2_i4 = inv.process(c2_i3);
  cout << "INV complete" << endl;
  SwiftImage<uint16> c2_i5 = wat.process(c2_i4);
  cout << "Wat complete" << endl;
  SwiftImage<uint16> c2_i6 = c2_i2 && c2_i5;
  vector<SwiftImageObject<> > c2_objs1 = segmenter.process(c2_i6);
  
 
  // Now have c1_objs1 and c2_objs1 for registration


  ChannelRegistration<> chanreg(c2_i6.image_width(),c2_i6.image_height());
 
  vector<vector<SwiftImageObject<> > > cycles;
  cycles.push_back(c1_objs1);
  cycles.push_back(c2_objs1);

  vector<vector<SwiftImageCluster<> > > clusters = chanreg.process_channel_registration(cycles,SwiftImageCluster<>::base_a);
  
  for(int n=0;n<clusters.size();n++) {
    ut.test(static_cast<int>(clusters[n].size()),2);
  }
  
  ut.test(static_cast<int>(clusters.size()),5);
  
  // Cluster 1, Cycle 1
  ut.test(static_cast<int>(clusters[0][0].features[SwiftImageCluster<>::base_a].size()),1);
  ut.test(static_cast<int>(clusters[0][0].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[0][0].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[0][0].features[SwiftImageCluster<>::base_c].size()),0);

  ut.test(clusters[0][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,1);
  ut.test(clusters[0][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,1);
  ut.test(clusters[0][0].features[SwiftImageCluster<>::base_a][0].pixels[0].length,2);
  
  ut.test(clusters[0][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,1);
  ut.test(clusters[0][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,2);
  ut.test(clusters[0][0].features[SwiftImageCluster<>::base_a][0].pixels[1].length,2);
  
  // Cluster 1, Cycle 2
  ut.test(static_cast<int>(clusters[0][1].features[SwiftImageCluster<>::base_a].size()),1);
  ut.test(static_cast<int>(clusters[0][1].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[0][1].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[0][1].features[SwiftImageCluster<>::base_c].size()),0);

  ut.test(clusters[0][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,1);
  ut.test(clusters[0][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,1);
  ut.test(clusters[0][1].features[SwiftImageCluster<>::base_a][0].pixels[0].length,2);
  
  ut.test(clusters[0][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,1);
  ut.test(clusters[0][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,2);
  ut.test(clusters[0][1].features[SwiftImageCluster<>::base_a][0].pixels[1].length,2);
  
  // Cluster 2, Cycle 1
  ut.test(static_cast<int>(clusters[1][0].features[SwiftImageCluster<>::base_a].size()),2);
  ut.test(static_cast<int>(clusters[1][0].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[1][0].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[1][0].features[SwiftImageCluster<>::base_c].size()),0);
  
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,1);
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,5);
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][0].pixels[0].length,2);
  
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,1);
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,6);
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][0].pixels[1].length,2);

  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][1].pixels[0].pos.x,5);
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][1].pixels[0].pos.y,5);
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][1].pixels[0].length,2);
  
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][1].pixels[1].pos.x,5);
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][1].pixels[1].pos.y,6);
  ut.test(clusters[1][0].features[SwiftImageCluster<>::base_a][1].pixels[1].length,2);

  // Cluster 2, Cycle 2
  ut.test(static_cast<int>(clusters[1][1].features[SwiftImageCluster<>::base_a].size()),1);
  ut.test(static_cast<int>(clusters[1][1].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[1][1].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[1][1].features[SwiftImageCluster<>::base_c].size()),0);

  ut.test(clusters[1][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,1);
  ut.test(clusters[1][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,5);
  ut.test(clusters[1][1].features[SwiftImageCluster<>::base_a][0].pixels[0].length,6);
  
  ut.test(clusters[1][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,1);
  ut.test(clusters[1][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,6);
  ut.test(clusters[1][1].features[SwiftImageCluster<>::base_a][0].pixels[1].length,6);
  
  // Cluster 3, Cycle 1
  ut.test(static_cast<int>(clusters[2][0].features[SwiftImageCluster<>::base_a].size()),1);
  ut.test(static_cast<int>(clusters[2][0].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[2][0].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[2][0].features[SwiftImageCluster<>::base_c].size()),0);
  
  ut.test(clusters[2][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,8);
  ut.test(clusters[2][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,5);
  ut.test(clusters[2][0].features[SwiftImageCluster<>::base_a][0].pixels[0].length,2);
  
  ut.test(clusters[2][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,8);
  ut.test(clusters[2][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,6);
  ut.test(clusters[2][0].features[SwiftImageCluster<>::base_a][0].pixels[1].length,2);
  
  // Cluster 3, Cycle 2
  ut.test(static_cast<int>(clusters[2][1].features[SwiftImageCluster<>::base_a].size()),1);
  ut.test(static_cast<int>(clusters[2][1].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[2][1].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[2][1].features[SwiftImageCluster<>::base_c].size()),0);
  
  ut.test(clusters[2][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,8);
  ut.test(clusters[2][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,5);
  ut.test(clusters[2][1].features[SwiftImageCluster<>::base_a][0].pixels[0].length,2);
  
  ut.test(clusters[2][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,8);
  ut.test(clusters[2][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,6);
  ut.test(clusters[2][1].features[SwiftImageCluster<>::base_a][0].pixels[1].length,2);
  
  // Cluster 4, Cycle 1
  ut.test(static_cast<int>(clusters[4][0].features[SwiftImageCluster<>::base_a].size()),1);
  ut.test(static_cast<int>(clusters[4][0].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[4][0].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[4][0].features[SwiftImageCluster<>::base_c].size()),0);

  ut.test(clusters[4][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,5); 
  ut.test(clusters[4][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,1);
  ut.test(clusters[4][0].features[SwiftImageCluster<>::base_a][0].pixels[0].length,2);
  
  ut.test(clusters[4][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,5);
  ut.test(clusters[4][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,2);
  ut.test(clusters[4][0].features[SwiftImageCluster<>::base_a][0].pixels[1].length,2);
 
  // Cluster 4, Cycle 2
  ut.test(static_cast<int>(clusters[4][1].features[SwiftImageCluster<>::base_a].size()),1);
  ut.test(static_cast<int>(clusters[4][1].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[4][1].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[4][1].features[SwiftImageCluster<>::base_c].size()),0);

  ut.test(clusters[4][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,5);
  ut.test(clusters[4][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,1);
  ut.test(clusters[4][1].features[SwiftImageCluster<>::base_a][0].pixels[0].length,2);
  
  ut.test(clusters[4][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,5);
  ut.test(clusters[4][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,2);
  ut.test(clusters[4][1].features[SwiftImageCluster<>::base_a][0].pixels[1].length,2);

  // Cluster 5, Cycle 1
  ut.test(static_cast<int>(clusters[3][0].features[SwiftImageCluster<>::base_a].size()),1);
  ut.test(static_cast<int>(clusters[3][0].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[3][0].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[3][0].features[SwiftImageCluster<>::base_c].size()),0);
  
  ut.test(clusters[3][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,1);
  ut.test(clusters[3][0].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,8);
  ut.test(clusters[3][0].features[SwiftImageCluster<>::base_a][0].pixels[0].length,9);
  
  ut.test(clusters[3][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,1);
  ut.test(clusters[3][0].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,9);
  ut.test(clusters[3][0].features[SwiftImageCluster<>::base_a][0].pixels[1].length,9);

 
  // Cluster 5, Cycle 2
  ut.test(static_cast<int>(clusters[3][1].features[SwiftImageCluster<>::base_a].size()),3);
  ut.test(static_cast<int>(clusters[3][1].features[SwiftImageCluster<>::base_t].size()),0);
  ut.test(static_cast<int>(clusters[3][1].features[SwiftImageCluster<>::base_g].size()),0);
  ut.test(static_cast<int>(clusters[3][1].features[SwiftImageCluster<>::base_c].size()),0);
  
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.x,1);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][0].pixels[0].pos.y,8);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][0].pixels[0].length,2);
  
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.x,1);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][0].pixels[1].pos.y,9);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][0].pixels[1].length,2);
  
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][1].pixels[0].pos.x,4);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][1].pixels[0].pos.y,8);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][1].pixels[0].length,2);

  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][1].pixels[1].pos.x,4);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][1].pixels[1].pos.y,9);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][1].pixels[1].length,2);
  
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][2].pixels[0].pos.x,7);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][2].pixels[0].pos.y,8);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][2].pixels[0].length,2);

  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][2].pixels[1].pos.x,7);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][2].pixels[1].pos.y,9);
  ut.test(clusters[3][1].features[SwiftImageCluster<>::base_a][2].pixels[1].length,2);
  
  ut.end_test_set();
}

