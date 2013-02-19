#ifndef UNITTEST_H
#define UNITTEST_H

#include <iostream>

#define test(t1,t2) _test(t1,t2,__LINE__,__FILE__)
#define test_approx(t1,t2,t3) _test_approx(t1,t2,t3,__LINE__,__FILE__)

class UnitTest {
private:
  int tests_failed;
  int tests_passed;
  int total_tests_failed;
  int total_tests_passed;
  std::string test_set_name;
  std::string current_description;

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
