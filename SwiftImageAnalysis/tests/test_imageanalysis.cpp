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
#include <vector>
#include "ImageAnalysis.h"
#include "Cluster.h"

using namespace std;


void test_registrationchannel(UnitTest &ut) {
  // Testing Cross-channel Registration
  vector<string> images_a_filenames;
  vector<string> images_c_filenames;
  vector<string> images_g_filenames;
  vector<string> images_t_filenames;  

  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C1.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C2.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C3.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C4.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C5.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C6.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C7.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C8.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake4/Images/L001/C9.1/s_1_1_a.tif");
  
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C1.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C2.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C3.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C4.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C5.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C6.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C7.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C8.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake4/Images/L001/C9.1/s_1_1_t.tif");
  
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C1.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C2.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C3.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C4.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C5.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C6.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C7.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C8.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake4/Images/L001/C9.1/s_1_1_g.tif");
  
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C1.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C2.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C3.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C4.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C5.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C6.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C7.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C8.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake4/Images/L001/C9.1/s_1_1_c.tif");

  ImageAnalysis<> m_image_analyser(images_a_filenames,images_c_filenames,images_g_filenames,images_t_filenames);
  m_image_analyser.params_threshold_window            = 6;
  m_image_analyser.params_threshold_correlation_window= 6;

  m_image_analyser.params_threshold             = 0.6;
  m_image_analyser.params_threshold_correlation = 0.6;
  
  m_image_analyser.params_agregate_cycle        = 1;
  m_image_analyser.params_reference_cycle       = 0;
  m_image_analyser.params_offsets_median        = false;
  m_image_analyser.params_combine_reference     = false;
  m_image_analyser.params_crosschannelsave      = true;
  m_image_analyser.params_correlation_method    = ChannelOffsets<uint16>::fft;

  vector<Cluster<> > clusters = m_image_analyser.generate();

  ut.test(static_cast<int>(clusters.size()),2);
  vector<ReadIntensity<> > &signal = clusters[0].signal("RAW");
 
  ut.test(signal[0].get_base(ReadIntensity<>::base_a),static_cast<double>(29.0));
  ut.test(signal[0].get_base(ReadIntensity<>::base_c),static_cast<double>(29.0));
  ut.test(signal[0].get_base(ReadIntensity<>::base_g),static_cast<double>(39.0));
  ut.test(signal[0].get_base(ReadIntensity<>::base_t),static_cast<double>(49.0));
  ut.test(signal[1].get_base(ReadIntensity<>::base_a),static_cast<double>(59.0));
  ut.test(signal[1].get_base(ReadIntensity<>::base_c),static_cast<double>(69.0));
  ut.test(signal[1].get_base(ReadIntensity<>::base_g),static_cast<double>(79.0));
  ut.test(signal[1].get_base(ReadIntensity<>::base_t),static_cast<double>(89.0));
  ut.test(signal[2].get_base(ReadIntensity<>::base_a),static_cast<double>(9.0));
  ut.test(signal[2].get_base(ReadIntensity<>::base_c),static_cast<double>(10.0));
  ut.test(signal[2].get_base(ReadIntensity<>::base_g),static_cast<double>(11.0));
  ut.test(signal[2].get_base(ReadIntensity<>::base_t),static_cast<double>(12.0));
  ut.test(signal[3].get_base(ReadIntensity<>::base_a),static_cast<double>(13.0));
  ut.test(signal[3].get_base(ReadIntensity<>::base_c),static_cast<double>(14.0));
  ut.test(signal[3].get_base(ReadIntensity<>::base_g),static_cast<double>(15.0));
  ut.test(signal[3].get_base(ReadIntensity<>::base_t),static_cast<double>(16.0));
  ut.test(signal[4].get_base(ReadIntensity<>::base_a),static_cast<double>(17.0));
  ut.test(signal[4].get_base(ReadIntensity<>::base_c),static_cast<double>(18.0));
  ut.test(signal[4].get_base(ReadIntensity<>::base_g),static_cast<double>(19.0));
  ut.test(signal[4].get_base(ReadIntensity<>::base_t),static_cast<double>(20.0));
  ut.test(signal[5].get_base(ReadIntensity<>::base_a),static_cast<double>(21.0));
  ut.test(signal[5].get_base(ReadIntensity<>::base_c),static_cast<double>(22.0));
  ut.test(signal[5].get_base(ReadIntensity<>::base_g),static_cast<double>(23.0));
  ut.test(signal[5].get_base(ReadIntensity<>::base_t),static_cast<double>(24.0));
  ut.test(signal[6].get_base(ReadIntensity<>::base_a),static_cast<double>(25.0));
  ut.test(signal[6].get_base(ReadIntensity<>::base_c),static_cast<double>(26.0));
  ut.test(signal[6].get_base(ReadIntensity<>::base_g),static_cast<double>(27.0));
  ut.test(signal[6].get_base(ReadIntensity<>::base_t),static_cast<double>(28.0));
  ut.test(signal[7].get_base(ReadIntensity<>::base_a),static_cast<double>(29.0));
  ut.test(signal[7].get_base(ReadIntensity<>::base_c),static_cast<double>(30.0));
  ut.test(signal[7].get_base(ReadIntensity<>::base_g),static_cast<double>(31.0));
  ut.test(signal[7].get_base(ReadIntensity<>::base_t),static_cast<double>(32.0));
  ut.test(signal[8].get_base(ReadIntensity<>::base_a),static_cast<double>(33.0));
  ut.test(signal[8].get_base(ReadIntensity<>::base_c),static_cast<double>(34.0));
  ut.test(signal[8].get_base(ReadIntensity<>::base_g),static_cast<double>(35.0));
  ut.test(signal[8].get_base(ReadIntensity<>::base_t),static_cast<double>(36.0));
}

void test_registrationcrosschannel(UnitTest &ut) {
  // Testing Cross-channel Registration
  vector<string> images_a_filenames;
  vector<string> images_c_filenames;
  vector<string> images_g_filenames;
  vector<string> images_t_filenames;  

  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C1.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C2.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C3.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C4.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C5.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C6.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C7.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C8.1/s_1_1_a.tif");
  images_a_filenames.push_back("./Images/smallfake3/Images/L001/C9.1/s_1_1_a.tif");
  
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C1.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C2.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C3.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C4.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C5.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C6.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C7.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C8.1/s_1_1_t.tif");
  images_t_filenames.push_back("./Images/smallfake3/Images/L001/C9.1/s_1_1_t.tif");
  
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C1.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C2.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C3.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C4.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C5.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C6.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C7.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C8.1/s_1_1_g.tif");
  images_g_filenames.push_back("./Images/smallfake3/Images/L001/C9.1/s_1_1_g.tif");
  
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C1.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C2.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C3.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C4.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C5.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C6.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C7.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C8.1/s_1_1_c.tif");
  images_c_filenames.push_back("./Images/smallfake3/Images/L001/C9.1/s_1_1_c.tif");


  ImageAnalysis<> m_image_analyser(images_a_filenames,images_c_filenames,images_g_filenames,images_t_filenames);
  m_image_analyser.params_threshold             = 0.6;
  m_image_analyser.params_threshold_correlation = 0.6;
  m_image_analyser.params_agregate_cycle        = 1;
  m_image_analyser.params_reference_cycle       = 1;
  m_image_analyser.params_offsets_median        = false;

  vector<Cluster<> > clusters = m_image_analyser.generate();

  vector<ReadIntensity<> > &signal = clusters[0].signal("RAW");

  ut.test(signal[0].get_base(ReadIntensity<>::base_a),static_cast<double>(1.0));
  ut.test(signal[0].get_base(ReadIntensity<>::base_c),static_cast<double>(2.0));
  ut.test(signal[0].get_base(ReadIntensity<>::base_g),static_cast<double>(3.0));
  ut.test(signal[0].get_base(ReadIntensity<>::base_t),static_cast<double>(4.0));
  ut.test(signal[1].get_base(ReadIntensity<>::base_a),static_cast<double>(5.0));
  ut.test(signal[1].get_base(ReadIntensity<>::base_c),static_cast<double>(6.0));
  ut.test(signal[1].get_base(ReadIntensity<>::base_g),static_cast<double>(7.0));
  ut.test(signal[1].get_base(ReadIntensity<>::base_t),static_cast<double>(8.0));
  ut.test(signal[2].get_base(ReadIntensity<>::base_a),static_cast<double>(9.0));
  ut.test(signal[2].get_base(ReadIntensity<>::base_c),static_cast<double>(10.0));
  ut.test(signal[2].get_base(ReadIntensity<>::base_g),static_cast<double>(11.0));
  ut.test(signal[2].get_base(ReadIntensity<>::base_t),static_cast<double>(12.0));
  ut.test(signal[3].get_base(ReadIntensity<>::base_a),static_cast<double>(13.0));
  ut.test(signal[3].get_base(ReadIntensity<>::base_c),static_cast<double>(14.0));
  ut.test(signal[3].get_base(ReadIntensity<>::base_g),static_cast<double>(15.0));
  ut.test(signal[3].get_base(ReadIntensity<>::base_t),static_cast<double>(16.0));
  ut.test(signal[4].get_base(ReadIntensity<>::base_a),static_cast<double>(17.0));
  ut.test(signal[4].get_base(ReadIntensity<>::base_c),static_cast<double>(18.0));
  ut.test(signal[4].get_base(ReadIntensity<>::base_g),static_cast<double>(19.0));
  ut.test(signal[4].get_base(ReadIntensity<>::base_t),static_cast<double>(20.0));
  ut.test(signal[5].get_base(ReadIntensity<>::base_a),static_cast<double>(21.0));
  ut.test(signal[5].get_base(ReadIntensity<>::base_c),static_cast<double>(22.0));
  ut.test(signal[5].get_base(ReadIntensity<>::base_g),static_cast<double>(23.0));
  ut.test(signal[5].get_base(ReadIntensity<>::base_t),static_cast<double>(24.0));
  ut.test(signal[6].get_base(ReadIntensity<>::base_a),static_cast<double>(25.0));
  ut.test(signal[6].get_base(ReadIntensity<>::base_c),static_cast<double>(26.0));
  ut.test(signal[6].get_base(ReadIntensity<>::base_g),static_cast<double>(27.0));
  ut.test(signal[6].get_base(ReadIntensity<>::base_t),static_cast<double>(28.0));
  ut.test(signal[7].get_base(ReadIntensity<>::base_a),static_cast<double>(29.0));
  ut.test(signal[7].get_base(ReadIntensity<>::base_c),static_cast<double>(30.0));
  ut.test(signal[7].get_base(ReadIntensity<>::base_g),static_cast<double>(31.0));
  ut.test(signal[7].get_base(ReadIntensity<>::base_t),static_cast<double>(32.0));
  ut.test(signal[8].get_base(ReadIntensity<>::base_a),static_cast<double>(33.0));
  ut.test(signal[8].get_base(ReadIntensity<>::base_c),static_cast<double>(34.0));
  ut.test(signal[8].get_base(ReadIntensity<>::base_g),static_cast<double>(35.0));
  ut.test(signal[8].get_base(ReadIntensity<>::base_t),static_cast<double>(36.0));
}

void test_imageanalysis(UnitTest &ut) {

  ut.begin_test_set("ImageAnalysis");

  test_registrationchannel(ut);
//  test_registrationcrosschannel(ut);
}
