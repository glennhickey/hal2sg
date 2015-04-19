/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include <sstream>
#include <ctime>
#include <cmath>
#include "unitTests.h"
#include "snphandler.h"

using namespace std;

void snpMapTest(CuTest *testCase)
{
}

CuSuite* snpHandlerTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, snpMapTest);
  return suite;
}
