#include "utf.h"
#include "test_intensity2filterblobs.h"
#include "intensity2filterblobs.h"

void test_intensity2filterblobs(UnitTest &ut) {

  ut.begin_test_set("intensity2filterblobs");

  intensity2filterblobs("testdata2.int","testdata2.nse","out_int","out_nse",10);

  // Very lazy, probably not very portable... honestly most of my tests are much better than this.
  int ret=-1;
  ret = system("diff out_int blobfilter.int.real");
  ut.test(ret,0);
  
  ret = system("diff out_nse blobfilter.nse.real");
  ut.test(ret,0);

  ut.end_test_set();
}

