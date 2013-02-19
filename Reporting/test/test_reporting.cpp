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
#include "Reporting.h"
#include "test_reporting.h"
#include "Cluster.h"
#include "ScoredSequence.h"

void test_reporting(UnitTest &ut) {

  ut.begin_test_set("Reporting");

  Cluster<double> c1;
  
  c1.add_signal("RAW");
  c1.signal("RAW").push_back(ReadIntensity<double>(1000 ,0,0,0   ));
  c1.signal("RAW").push_back(ReadIntensity<double>(500  ,0,0,2000));
  c1.signal("RAW").push_back(ReadIntensity<double>(250  ,0,0,2000));
  c1.signal("RAW").push_back(ReadIntensity<double>(125  ,0,0,2000));
  c1.signal("RAW").push_back(ReadIntensity<double>(62.5 ,0,0,2000));
  c1.signal("RAW").push_back(ReadIntensity<double>(31.25,0,0,2000));
  
  c1.add_signal("FINAL_CORRECTED");
  c1.signal("FINAL_CORRECTED").push_back(ReadIntensity<double>(1000 ,0,0,0   ));
  c1.signal("FINAL_CORRECTED").push_back(ReadIntensity<double>(500  ,0,0,2000));
  c1.signal("FINAL_CORRECTED").push_back(ReadIntensity<double>(250  ,0,0,2000));
  c1.signal("FINAL_CORRECTED").push_back(ReadIntensity<double>(125  ,0,0,2000));
  c1.signal("FINAL_CORRECTED").push_back(ReadIntensity<double>(62.5 ,0,0,2000));
  c1.signal("FINAL_CORRECTED").push_back(ReadIntensity<double>(31.25,0,0,2000));

  c1.add_sequence("BASECALL_FINAL");
  //TODO: I'm really unclear as to why I need this cast, for some reason gcc can't see base_type from here...
  c1.sequence("BASECALL_FINAL").sequence().push_back(static_cast<int>(ScoredSequence<>::base_a));
  c1.sequence("BASECALL_FINAL").sequence().push_back(static_cast<int>(ScoredSequence<>::base_t));
  c1.sequence("BASECALL_FINAL").sequence().push_back(static_cast<int>(ScoredSequence<>::base_t));
  c1.sequence("BASECALL_FINAL").sequence().push_back(static_cast<int>(ScoredSequence<>::base_t));
  c1.sequence("BASECALL_FINAL").sequence().push_back(static_cast<int>(ScoredSequence<>::base_t));
  c1.sequence("BASECALL_FINAL").sequence().push_back(static_cast<int>(ScoredSequence<>::base_t));
  
  vector<Cluster<double> > clusters;
  clusters.push_back(c1);
  clusters.push_back(c1);
  clusters.push_back(c1);
  clusters.push_back(c1);
  clusters.push_back(c1);
  clusters.push_back(c1);

  Reporting<double> m_reporting(clusters,true,"./phi_plus_SNPs.fa");

  m_reporting.write_report_file("report");

  ut.end_test_set();
}

