#include "utf.h"
#include "test_intensity2pairwiseplot.h"
#include "test_readintensity.h"
#include "test_intensitysequence.h"
#include "test_intensity2filterblobs.h"

int main(void) {

  UnitTest ut("Intensity Tools");

  test_intensity2pairwiseplot(ut);  
  test_readintensity(ut);  
  test_intensitysequence(ut);  
  test_intensity2filterblobs(ut);

  ut.test_report();

  return 0;
}

