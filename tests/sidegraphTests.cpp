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
  SGPosition p1;
  SGPosition p2;
  SGPosition p3;

  p1._seqid = 0;
  p1._pos = 1;
  p2._seqid = 1;
  p2._pos = 0;
  p3._seqid = 0;
  p3._pos = 0;
  CuAssertTrue(testCase, p1 == p1);
  CuAssertTrue(testCase, p1 < p2);
  CuAssertTrue(testCase, p3 < p1);

  SGSide s1;
  SGSide s2;
  SGSide s3;
  s1._base = p1;
  s1._forward = true;
  s2._base = p2;
  s2._forward = false;
  s3._base = p1;
  s3._forward = false;
  CuAssertTrue(testCase, s1 == s1);
  CuAssertTrue(testCase, s1 < s2);
  CuAssertTrue(testCase, s3 < s1);

  SGJoin j1;
  SGJoin j2;
  SGJoin j3;
  j1._side1 = s1;
  j1._side2 = s2;
  j2._side1 = s2;
  j2._side2 = s3;
  j3._side1 = s1;
  j3._side2 = s3;
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
    s->_length = i * 10;
    stringstream ss;
    ss << "seq_" << i; 
    s->_name = ss.str();
    sg.addSequence(s);
  }
  for (size_t i = 0; i < 1000; ++i)
  {
    const SGSequence* s = sg.getSequence(i);
    CuAssertTrue(testCase, s->_id == i);
    CuAssertTrue(testCase, s->_length == i * 10);
    stringstream ss;
    ss << "seq_" << i;
    CuAssertTrue(testCase, s->_name == ss.str());
  }
  
  for (size_t i = 0; i < 100; ++i)
  {
    for (size_t j = 0; j < 100; ++j)
    {
      SGJoin* join = new SGJoin();
      join->_side1._forward = true;
      join->_side1._base._pos = i;
      join->_side1._base._seqid = i;
      join->_side2._forward = false;
      join->_side2._base._pos = j + 100;
      join->_side2._base._seqid = i;
      CuAssertTrue(testCase, sg.addJoin(join) == join);      
    }
  }
  for (size_t i = 0; i < 100; ++i)
  {
    for (size_t j = 0; j < 100; ++j)
    {
      SGJoin* join = new SGJoin();
      join->_side1._forward = true;
      join->_side1._base._pos = i;
      join->_side1._base._seqid = i;
      join->_side2._forward = false;
      join->_side2._base._pos = j + 100;
      join->_side2._base._seqid = i;
      const SGJoin* join2 = sg.getJoin(join);      
      CuAssertTrue(testCase, *join == *join2);
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
