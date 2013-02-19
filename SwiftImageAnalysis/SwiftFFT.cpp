/*
    Swift (c) 2008 Genome Research Ltd.
    Authors: Nava Whiteford and Tom Skelly (new@sgenomics.org ts6@sanger.ac.uk)

    This file is part of Swift (http://swiftng.sourceforge.net).

    Swift is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Swift is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Swift.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "SwiftFFT.h"

// What FFT batch size should we choose? According to the user guide, "FFTW is
// best at handling sizes of the form 2^a*3^b*5^c*7^d*11^e*13^f, where e+f is
// either 0 or 1, and the other exponents are arbitrary. Other sizes are computed
// by means of a slow, general-purpose algorithm (which nevertheless retains 
// O(n log n)."
//
// We'll keep it simple and use only sizes which are powers of 2 and/or 3. The table below 
// contains all of the possibilities. In order for the DC-in-the-middle hack to work,
// only even sizes are present -- no pure powers of 3. We'll zero-pad the input
// data up to the required size. [NOTE: the DC-in-the-middle hack has been disabled --
// couldn't get it working, and it was only going to be useful for debugging.
//
// The wisdom of that is borne out in practice: A 1004*1002 FFT takes 16 seconds to plan,
// and 200 msec to execute. A 1024*1024 FFT takes around a second to plan, and 60 msec to
// execute.

const int SwiftFFT::good_sizes[] = {
      512,      576,      648,      768,      864,      972,     1024,     1152,
     1296,     1458,     1536,     1728,     1944,     2048,     2304,     2592,
     2916,     3072,     3456,     3888,     4096,     4374,     4608,     5184,
     5832,     6144,     6912,     7776,     8192,     8748,     9216,    10368,
    11664,    12288,    13122,    13824,    15552,    16384,    17496,    18432,
    20736,    23328,    24576,    26244,    27648,    31104,    32768,    34992,
    36864,    39366,    41472,    46656,    49152,    52488,    55296,    62208,
    65536,    69984,    73728,    78732,    82944,    93312,    98304,   104976,
   110592,   118098,   124416,   131072,   139968,   147456,   157464,   165888,
   186624,   196608,   209952,   221184,   236196,   248832,   262144,   279936,
   294912,   314928,   331776,   354294,   373248,   393216,   419904,   442368,
   472392,   497664,   524288,   559872,   589824,   629856,   663552,   708588,
   746496,   786432,   839808,   884736,   944784,   995328,  1048576,  1062882,
  1119744,  1179648,  1259712,  1327104,  1417176,  1492992,  1572864,  1679616,
  1769472,  1889568,  1990656,  2097152,  2125764,  2239488,  2359296,  2519424,
  2654208,  2834352,  2985984,  3145728,  3188646,  3359232,  3538944,  3779136,
  3981312,  4194304,  4251528,  4478976,  4718592,  5038848,  5308416,  5668704,
  5971968,  6291456,  6377292,  6718464,  7077888,  7558272,  7962624,  8388608,
  8503056,  8957952,  9437184,  9565938, 10077696, 10616832, 11337408, 11943936,
 12582912, 12754584, 13436928, 14155776, 15116544, 15925248, 16777216
};

const int SwiftFFT::num_sizes = sizeof(good_sizes) / sizeof(int);

int SwiftFFT::pick_size (int size) {
  int i(0);
  while (i<num_sizes && good_sizes[i] < size) {
    ++i;
  }
  if (i >= num_sizes) {
    throw(std::out_of_range("FFT size too large"));
  }
  return (good_sizes[i]);
}
  
