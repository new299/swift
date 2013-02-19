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

#ifndef UNITTEST_H
#define UNITTEST_H

#include <iostream>

#define test(t1,t2) _test(t1,t2,__LINE__,__FILE__)
#define test_gt(t1,t2) _test_gt(t1,t2,__LINE__,__FILE__)
#define test_equality(t1,t2) _test(t1,t2,__LINE__,__FILE__)
#define test_truth(t1) _test(t1,__LINE__,__FILE__)
#define test_approx(t1,t2,t3) _test_approx(t1,t2,t3,__LINE__,__FILE__)

class UnitTest {
private:
  int tests_failed;
  int tests_passed;
  int total_tests_failed;
  int total_tests_passed;
  std::string current_description;
  std::string test_set_name;

public:

  UnitTest(std::string test_set_name_in) : tests_failed(0),
                                           tests_passed(0),
                                           total_tests_failed(0),
                                           total_tests_passed(0),
                                           current_description(),
                                           test_set_name(test_set_name_in) {
    std::cout << "*** Test set  : " << test_set_name << std::endl;
  }

  void begin_test_set(std::string description) {
    current_description = description;
    tests_failed = 0;
    tests_passed = 0;
    std::cout << "****** Testing: " << current_description << std::endl;
  }

  void end_test_set() {
    std::cout << "****** Test   : " << current_description << " complete, ";
    std::cout << "passed " << tests_passed;
    std::cout << ", failed " << tests_failed << "." << std::endl;
  }

  void _test(bool t1,int linenumber,const char *current_file) {
    if(!t1) {
      std::cout << "****** FAILED : " << current_file << "," << linenumber << std::endl;
      total_tests_failed++;
      tests_failed++;
    } else { tests_passed++; total_tests_passed++; }
  }

  template<class _TestType>
  void _test(_TestType t1,_TestType t2,int linenumber,const char *current_file) {
    bool test_result = (t1 == t2);

    if(!test_result) {
        std::cout << "****** FAILED : " << current_file << "," << linenumber;
	std::cout << ": " << t1 << " is not equal to " << t2 << std::endl;
        total_tests_failed++;
        tests_failed++;
    } else { tests_passed++; total_tests_passed++; }
  }
  
  template<class _TestType>
  void _test_gt(_TestType t1,_TestType t2,int linenumber,const char *current_file) {
    bool test_result = (t1 > t2);

    if(!test_result) {
        std::cout << "****** FAILED : " << current_file << "," << linenumber;
	std::cout << ": " << t1 << " is not greater than " << t2 << std::endl;
        total_tests_failed++;
        tests_failed++;
    } else { tests_passed++; total_tests_passed++; }
  }

  template<class _TestType>
  void _test_approx(_TestType t1,_TestType t2,_TestType max_error,int linenumber,const char *current_file) {
    _TestType difference = t1-t2;
    if(difference < 0) difference = 0-difference;

    if(difference > max_error) {
      std::cout << "****** FAILED : " << current_file << "," << linenumber;
      std::cout << ": " << t1 << " is not equal to " << t2;
      std::cout << " within difference of " << max_error << std::endl;
      tests_failed++;
      total_tests_failed++;
    } else { tests_passed++; total_tests_passed++; }
  }

  void test_report() {
    std::cout << "*** Test set  : " << test_set_name << " complete, ";
    std::cout << "passed " << total_tests_passed;
    std::cout << " failed " << total_tests_failed << "." << std::endl;
    if(total_tests_failed != 0) std::cout << "*** TEST FAILED!" << std::endl;
  }

};

#endif
