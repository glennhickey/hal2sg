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

/** test lengthTo function */
void sideGraphTestSide(CuTest *testCase)
{
  SGPosition p1(0, 5);
  SGPosition p2(0, 10);
  SGPosition p3(1, 1);

  SGSide s1f(p1, false);
  SGSide s1r(p1, true);
  SGSide s2f(p2, false);
  SGSide s2r(p2, true);
  SGSide s3f(p3, false);
  SGSide s3r(p3, true);

  CuAssertTrue(testCase, s1f.lengthTo(s1f) == 0);
  CuAssertTrue(testCase, s1f.lengthTo(s1r) == 1);
  CuAssertTrue(testCase, s1f.lengthTo(s2f) == 5);
  CuAssertTrue(testCase, s1f.lengthTo(s2r) == 4);
  CuAssertTrue(testCase, s1f.lengthTo(s3f) == -1);
  CuAssertTrue(testCase, s1f.lengthTo(s3r) == -1);

  CuAssertTrue(testCase, s1r.lengthTo(s1f) == 1);
  CuAssertTrue(testCase, s1r.lengthTo(s1r) == 0);
  CuAssertTrue(testCase, s1r.lengthTo(s2f) == 6);
  CuAssertTrue(testCase, s1r.lengthTo(s2r) == 5);
  CuAssertTrue(testCase, s1r.lengthTo(s3f) == -1);
  CuAssertTrue(testCase, s1r.lengthTo(s3r) == -1);

  CuAssertTrue(testCase, s2f.lengthTo(s1f) == 5);
  CuAssertTrue(testCase, s2f.lengthTo(s1r) == 6);
  CuAssertTrue(testCase, s2f.lengthTo(s2f) == 0);
  CuAssertTrue(testCase, s2f.lengthTo(s2r) == 1);
  CuAssertTrue(testCase, s2f.lengthTo(s3f) == -1);
  CuAssertTrue(testCase, s2f.lengthTo(s3r) == -1);
  
  CuAssertTrue(testCase, s2r.lengthTo(s1f) == 4);
  CuAssertTrue(testCase, s2r.lengthTo(s1r) == 5);
  CuAssertTrue(testCase, s2r.lengthTo(s2f) == 1);
  CuAssertTrue(testCase, s2r.lengthTo(s2r) == 0);
  CuAssertTrue(testCase, s2r.lengthTo(s3f) == -1);
  CuAssertTrue(testCase, s2r.lengthTo(s3r) == -1);
}

/** test Segment construction from sides function */
void sideGraphTestSegment(CuTest *testCase)
{
  SGPosition p1(0, 5);
  SGPosition p2(0, 10);
  SGPosition p3(1, 1);

  SGSide s1f(p1, false);
  SGSide s1r(p1, true);
  SGSide s2f(p2, false);
  SGSide s2r(p2, true);
  SGSide s3f(p3, false);
  SGSide s3r(p3, true);

  SGSegment seg;
  SGPosition pos;
  
  seg = SGSegment(s1f, s1f);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s1f.getBase().getPos() + 1);
  CuAssertTrue(testCase, seg.getLength() == 0);

  seg = SGSegment(s1r, s1r);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s1r.getBase().getPos() - 1);
  CuAssertTrue(testCase, seg.getLength() == 0);

  seg = SGSegment(s1r, s1f);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s1r.getBase().getPos());
  CuAssertTrue(testCase, seg.getLength() == 1);

  seg = SGSegment(s1f, s1r);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s1r.getBase().getPos());
  CuAssertTrue(testCase, seg.getLength() == 1);

  seg = SGSegment(s1f, s2f);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s1f.getBase().getPos() + 1);
  CuAssertTrue(testCase, seg.getLength() == 5);

  seg = SGSegment(s1f, s2r);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s1f.getBase().getPos() + 1);
  CuAssertTrue(testCase, seg.getLength() == 4);

  seg = SGSegment(s1r, s2f);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s1f.getBase().getPos());
  CuAssertTrue(testCase, seg.getLength() == 6);

  seg = SGSegment(s1r, s2r);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s1f.getBase().getPos());
  CuAssertTrue(testCase, seg.getLength() == 5);

  seg = SGSegment(s2f, s1f);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s2f.getBase().getPos());
  CuAssertTrue(testCase, seg.getLength() == 5);

  seg = SGSegment(s2f, s1r);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s2f.getBase().getPos());
  CuAssertTrue(testCase, seg.getLength() == 6);

  seg = SGSegment(s2r, s1f);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s2f.getBase().getPos() - 1);
  CuAssertTrue(testCase, seg.getLength() == 4);

  seg = SGSegment(s2r, s1r);
  pos = seg.getSide().getBase();
  CuAssertTrue(testCase, pos.getPos() == s2f.getBase().getPos() - 1);
  CuAssertTrue(testCase, seg.getLength() == 5);
}


/** make sure basic structs sort okay */
void sideGraphTestJoin(CuTest *testCase)
{
  SGPosition p1(0, 1);
  SGPosition p2(1, 0);
  SGPosition p3(0, 0);

  CuAssertTrue(testCase, p1 == p1);
  CuAssertTrue(testCase, p1 < p2);
  CuAssertTrue(testCase, p3 < p1);

  SGSide s1(p1, false);
  SGSide s2(p2, true);
  SGSide s3(p1, true);

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
      SGJoin* join = new SGJoin(SGSide(SGPosition(i, i), false),
                                SGSide(SGPosition(j + 100, i), true));
      CuAssertTrue(testCase, sg.addJoin(join) == join);      
    }
  }
  for (size_t i = 0; i < 100; ++i)
  {
    for (size_t j = 0; j < 100; ++j)
    {
      SGJoin join(SGSide(SGPosition(i, i), false),
                  SGSide(SGPosition(j + 100, i), true));
      const SGJoin* join2 = sg.getJoin(&join);      
      CuAssertTrue(testCase, join == *join2);
    }
  }
}

CuSuite* sideGraphTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, sideGraphTestSide);
  SUITE_ADD_TEST(suite, sideGraphTestSegment);
  SUITE_ADD_TEST(suite, sideGraphTestJoin);
  SUITE_ADD_TEST(suite, sideGraphTest);
  return suite;
}
