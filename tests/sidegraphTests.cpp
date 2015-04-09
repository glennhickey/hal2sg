/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include <sstream>
#include "unitTests.h"
#include "sidegraph.h"

using namespace std;

/** make sure basic structs sort okay */
void sideGraphTestJoin(CuTest *testCase)
{
  SGPosition p1(0, 1);
  SGPosition p2(1, 0);
  SGPosition p3(0, 0);

  CuAssertTrue(testCase, p1 == p1);
  CuAssertTrue(testCase, p1 < p2);
  CuAssertTrue(testCase, p3 < p1);

  SGSide s1(p1, true);
  SGSide s2(p2, false);
  SGSide s3(p1, false);

  CuAssertTrue(testCase, s1 == s1);
  CuAssertTrue(testCase, s1 < s2);
  CuAssertTrue(testCase, s3 < s1);

  SGJoin j1(s1, s2);
  SGJoin j2(s2, s3);
  SGJoin j3(s1, s3);

  CuAssertTrue(testCase, j1 == j1);
  CuAssertTrue(testCase, j1 < j2);
  CuAssertTrue(testCase, j3 < j1);
  CuAssertTrue(testCase, SGJoinPtrLess().operator()(&j3, &j1));
}

/** basic get and sets for side graph */
void sideGraphTest(CuTest *testCase)
{
  SideGraph sg;
  // add some nonsense
  for (size_t i = 0; i < 1000; ++i)
  {
    SGSequence* s = new SGSequence();
    s->setLength(i * 10);
    stringstream ss;
    ss << "seq_" << i; 
    s->setName(ss.str());
    sg.addSequence(s);
  }
  for (size_t i = 0; i < 1000; ++i)
  {
    const SGSequence* s = sg.getSequence(i);
    CuAssertTrue(testCase, s->getID() == i);
    CuAssertTrue(testCase, s->getLength() == i * 10);
    stringstream ss;
    ss << "seq_" << i;
    CuAssertTrue(testCase, s->getName() == ss.str());
  }
  
  for (size_t i = 0; i < 100; ++i)
  {
    for (size_t j = 0; j < 100; ++j)
    {
      SGJoin* join = new SGJoin(SGSide(SGPosition(i, i), true),
                                SGSide(SGPosition(j + 100, i), false));
      CuAssertTrue(testCase, sg.addJoin(join) == join);      
    }
  }
  for (size_t i = 0; i < 100; ++i)
  {
    for (size_t j = 0; j < 100; ++j)
    {
      SGJoin join(SGSide(SGPosition(i, i), true),
                  SGSide(SGPosition(j + 100, i), false));
      const SGJoin* join2 = sg.getJoin(&join);      
      CuAssertTrue(testCase, join == *join2);
    }
  }
}

CuSuite* sideGraphTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, sideGraphTestJoin);
  SUITE_ADD_TEST(suite, sideGraphTest);
  return suite;
}
