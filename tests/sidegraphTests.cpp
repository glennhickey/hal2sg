/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include "unitTests.h"

void sideGraphTestSegments(CuTest *testCase)
{
}

CuSuite* sideGraphTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, sideGraphTestSegments);
  return suite;
}
