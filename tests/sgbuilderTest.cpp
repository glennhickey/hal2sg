/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include <sstream>
#include <ctime>
#include <cmath>
#include "halAlignmentTest.h"
#include "unitTests.h"
#include "sgbuilder.h"

using namespace std;
using namespace hal;

struct halTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

void halTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  Genome* leaf1Genome = alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);

  vector<Sequence::Info> seqVec(1);
  seqVec[0] =Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);
}

void halTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("AncGenome");
  const MetaData* ancMeta = ancGenome->getMetaData();
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
  CuAssertTrue(_testCase, leaf1Genome->getName() == "Leaf1");
  CuAssertTrue(_testCase, ancGenome->getSequenceLength() == 1000000);
  CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == 5000);
}

void sgBuilderHalTest(CuTest *testCase)
{
  try
  {
    halTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}
  
CuSuite* sgBuildTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, sgBuilderHalTest);
  return suite;
}
