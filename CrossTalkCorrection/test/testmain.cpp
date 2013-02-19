#include "utf.h"
#include "test_crosstalkcorrection.h"

int main(void) {

  UnitTest ut("CrossTalkCorrection Classes");

  test_crosstalkcorrection(ut);  
  
  ut.test_report();

  return 0;
}

