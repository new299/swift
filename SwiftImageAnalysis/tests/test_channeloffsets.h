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

#ifndef TEST_CHANNELOFFSETS
#define TEST_CHANNELOFFSETS

#include "ChannelOffsets.h"
#include "utf.h"
#include "SwiftImage.h"
#include "SwiftImageCluster.h"
#include "SwiftImageObject.h"
#include "Segmentation.h"
#include "CrossChannelRegistration.h"
#include "ChannelRegistration.h"
#include "RunLengthEncode.h"
#include "EuclideanDistanceMap.h"
#include "RLERun.h"
#include "Watershed.h"
#include "NWThreshold.h"
#include "Invert.h"
#include <math.h>

void test_channeloffsets(UnitTest &ut); 

#endif
