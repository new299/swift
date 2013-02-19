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
#include "ReadIntensity.h"
#include "test_readintensity.h"
#include "Cluster.h"

void test_readintensity(UnitTest &ut) {

  ut.begin_test_set("ReadIntensity");

  ReadIntensity<double> r(0,0,0,1000);

  ut.test(r.purity(),static_cast<double>(1));
 
  ReadIntensity<double> r1(0,0,500,1000);
  ut.test_approx(r1.purity(),0.66666666,0.0001);

  ReadIntensity<double> r2(0,10,500,1000);
  ut.test_approx(r2.purity(),0.66666666,0.0001);

  ReadIntensity<double> r3(0,600,1000,500);
  ut.test_approx(r3.purity(),0.625,0.0001);

  ReadIntensity<double> r4(1000,600,10,500);
  ut.test_approx(r4.purity(),0.625,0.0001);


  ut.test(r.max_base() ,ReadIntensity<double>::base_t);
  ut.test(r1.max_base(),ReadIntensity<double>::base_t);
  ut.test(r2.max_base(),ReadIntensity<double>::base_t);
  ut.test(r3.max_base(),ReadIntensity<double>::base_g);
  ut.test(r4.max_base(),ReadIntensity<double>::base_a);

  ut.test(r .sub_max_base() ,ReadIntensity<double>::base_g);
  ut.test(r1.sub_max_base() ,ReadIntensity<double>::base_g);
  ut.test(r2.sub_max_base() ,ReadIntensity<double>::base_g);
  ut.test(r3.sub_max_base() ,ReadIntensity<double>::base_c);
  ut.test(r4.sub_max_base() ,ReadIntensity<double>::base_c);

  ut.end_test_set();
}

