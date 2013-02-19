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

#include "utf.h"
#include "PhasingCorrection.h"
#include "test_phasingcorrection.h"
#include "Cluster.h"

void test_phasingcorrection(UnitTest &ut) {

  ut.begin_test_set("PhasingCorrection");

  Cluster<double> c1;
  c1.add_signal("RAW");
  c1.signal("RAW").push_back(ReadIntensity<double>(1000 ,0,0,0   ));
  c1.signal("RAW").push_back(ReadIntensity<double>(500  ,0,0,2000));
  c1.signal("RAW").push_back(ReadIntensity<double>(250  ,0,0,2000));
  c1.signal("RAW").push_back(ReadIntensity<double>(125  ,0,0,2000));
  c1.signal("RAW").push_back(ReadIntensity<double>(62.5 ,0,0,2000));
  c1.signal("RAW").push_back(ReadIntensity<double>(31.25,0,0,2000));

  vector<Cluster<double> > clusters;
  clusters.push_back(c1);

  PhasingCorrection<double> m_phasing_correction(0,0.5,0.6,"RAW","PHASE_CORRECTED");
  m_phasing_correction.process(clusters);
  
  ut.test(clusters[0].signal("PHASE_CORRECTED")[0],ReadIntensity<double>(1968.75,0,0,0));
  ut.test(clusters[0].signal("PHASE_CORRECTED")[1],ReadIntensity<double>(0      ,0,0,2000));
  ut.test(clusters[0].signal("PHASE_CORRECTED")[2],ReadIntensity<double>(0      ,0,0,2000));
  ut.test(clusters[0].signal("PHASE_CORRECTED")[3],ReadIntensity<double>(0      ,0,0,2000));
  ut.test(clusters[0].signal("PHASE_CORRECTED")[4],ReadIntensity<double>(0      ,0,0,2000));
  ut.test(clusters[0].signal("PHASE_CORRECTED")[5],ReadIntensity<double>(0      ,0,0,2000));

  ut.end_test_set();
}

