#include "utf.h"
#include "test_intensitysequence.h"
#include "intensity.h"

void test_intensitysequence(UnitTest &ut) {

  ut.begin_test_set("IntensitySequence");

  IntensitySequence is;

  ReadIntensity i0(1,2,3,4);
  ReadIntensity i1(1,2,3,4);
  ReadIntensity i2(1,2,3,4);
  ReadIntensity i3(1,2,3,4);
  ReadIntensity i4(1,2,3,4);
  ReadIntensity i5(1,2,3,4);
  ReadIntensity i6(1,2,3,4);
  ReadIntensity i7(1,2,3,4);
  ReadIntensity i8(1,2,3,4);
  ReadIntensity i9(1,2,3,4);
  ReadIntensity i10(1,2,3,4);
  ReadIntensity i11(1,2,3,4);
  ReadIntensity i12(1,2,3,4.2);
  ReadIntensity i13(1,2,3,3.2);
  ReadIntensity i14(1,2,3,3);
  ReadIntensity i15(1,2,3,3.1);

  is.intensities.push_back(i0);
  is.intensities.push_back(i1);
  is.intensities.push_back(i2);
  is.intensities.push_back(i3);
  is.intensities.push_back(i4);
  is.intensities.push_back(i5);
  is.intensities.push_back(i6);
  is.intensities.push_back(i7);
  is.intensities.push_back(i8);
  is.intensities.push_back(i9);
  is.intensities.push_back(i10);
  is.intensities.push_back(i11);
  is.intensities.push_back(i12);
  is.intensities.push_back(i13);
  is.intensities.push_back(i14);
  is.intensities.push_back(i15);

  ut.test(is.min_purity(1,11),(static_cast<double>(4)/static_cast<double>(4+3)));
  ut.test(is.min_purity(1,12),(static_cast<double>(4)/static_cast<double>(4+3)));
  ut.test(is.min_purity(1,13),(static_cast<double>(3.2)/static_cast<double>(3.2+3)));
  ut.test(is.min_purity(1,14),(static_cast<double>(3)/static_cast<double>(3+3)));

  ut.end_test_set();
}

