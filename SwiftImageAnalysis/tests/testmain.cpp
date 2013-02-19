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

#include "test_swiftimage.h"
#include "test_lowercomplete.h"
#include "test_morphologicalopening.h"
#include "test_morphologicalclosing.h"
#include "test_sobeloperator.h"
#include "test_adaptivethreshold.h"
#include "test_nwthreshold.h"
#include "test_euclideandistancemap.h"
#include "test_localmaxima.h"
#include "test_watershed.h"
#include "test_runlengthencode.h"
#include "utf.h"
#include "test_segmentation.h"
#include "test_channelregistration.h"
#include "test_crosschannelregistration.h"
#include "test_imageanalysis.h"
#include "test_channeloffsets.h"

int main(void) {

  UnitTest ut("Swift Image Analysis Classes");

  // test_lowercomplete(ut);
  test_channeloffsets(ut);
  //test_runlengthencode(ut);
  // test_segmentation(ut);
  //test_swiftimage(ut);  
  //test_crosschannelregistration(ut);
  //test_imageanalysis(ut);
  // test_channelregistration(ut);
  // test_watershed(ut);

  //test_morphologicalopening(ut);  
  //test_morphologicalclosing(ut);  
  //test_sobeloperator(ut);
  //test_adaptivethreshold(ut);
  //test_nwthreshold(ut);
  // test_euclideandistancemap(ut);
  // test_localmaxima(ut);

  ut.test_report();

  return 0;
}

