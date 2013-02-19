#include "utf.h"
#include "test_readintensity.h"
#include "intensity.h"

void test_readintensity(UnitTest &ut) {

  ut.begin_test_set("ReadIntensity");

  ReadIntensity i(1,2,3,4);

  ut.test(i.max_base(),ReadIntensity::base_c);
  ut.test(i.sub_max_base(),ReadIntensity::base_g);
  ut.test(i.getbase(i.max_base()),static_cast<double>(4));
  ut.test(i.getbase(i.sub_max_base()),static_cast<double>(3));
  ut.test(i.purity(),(static_cast<double>(4)/static_cast<double>(4+3)));
  
  ReadIntensity i1(5.2,34,3,2);

  ut.test(i1.max_base(),ReadIntensity::base_t);
  ut.test(i1.sub_max_base(),ReadIntensity::base_a);
  ut.test(i1.getbase(i1.max_base()),static_cast<double>(34));
  ut.test(i1.getbase(i1.sub_max_base()),static_cast<double>(5.2));
  ut.test(i1.purity(),(static_cast<double>(34)/static_cast<double>(34+5.2)));

  ut.end_test_set();
}

