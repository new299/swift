#include "utf.h"
#include "test_intensity2pairwiseplot.h"
#include "intensity2pairwiseplot.h"

void test_intensity2pairwiseplot(UnitTest &ut) {

  ut.begin_test_set("intensity2pairwiseplot");

  intensity2pairwiseplot("testdata","plot_");

  // Very lazy, probably not very portable... honestly most of my tests are much better than this.
  int ret=-1;
  ret = system("diff plot_ca plot_ca.real");
  ut.test(ret,0);
  
  ret = system("diff plot_cg plot_cg.real");
  ut.test(ret,0);
  
  ret = system("diff plot_ct plot_ct.real");
  ut.test(ret,0);
  
  ret = system("diff plot_ga plot_ga.real");
  ut.test(ret,0);
  
  ret = system("diff plot_gt plot_gt.real");
  ut.test(ret,0);

  ut.end_test_set();
}

