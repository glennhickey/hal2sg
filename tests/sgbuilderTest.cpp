/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <ctime>
#include <cmath>
#include "halAlignmentTest.h"
#include "unitTests.h"
#include "sgbuilder.h"

using namespace std;
using namespace hal;

///////////////////////////////////////////////////////////////////////////
//
//           BASIC HAL API INSTALLATION TEST
//
///////////////////////////////////////////////////////////////////////////

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
  alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);

  vector<Sequence::Info> seqVec(1);
  seqVec[0] =Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);
}

void halTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("AncGenome");
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

///////////////////////////////////////////////////////////////////////////
//
//           SINGLE INVERSION TEST
//
///////////////////////////////////////////////////////////////////////////

struct InversionTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

static string randDNA(size_t len)
{
  string out;
  for (size_t i = 0; i < len; ++i)
  {
    long x = rand() % 4;
    char c;
    switch(x)
    {
    case 0 : c = 'a'; break;
    case 1 : c = 'c'; break;
    case 2 : c = 'g'; break;
    default:
      c = 't';
    }
    out += c;
  }
  return out;
}

void InversionTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  Genome* leaf1Genome = alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);

  vector<Sequence::Info> seqVec(1);
  seqVec[0] =Sequence::Info("AncSequence", 100, 0, 10);
  ancGenome->setDimensions(seqVec);

  seqVec[0] =Sequence::Info("LeafSequence", 100, 10, 0);
  leaf1Genome->setDimensions(seqVec);

  string dna = randDNA(ancGenome->getSequenceLength());
  ancGenome->setString(dna);
  string leafdna = dna.substr(0, 50);
  string temp = dna.substr(50, 10);
  reverseComplement(temp);
  leafdna += temp;
  leafdna += dna.substr(60, 40);
  leaf1Genome->setString(leafdna);

  TopSegmentIteratorPtr top = leaf1Genome->getTopSegmentIterator();
  BottomSegmentIteratorPtr bottom = ancGenome->getBottomSegmentIterator();

  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    bottom->setTopParseIndex(NULL_INDEX);
    bottom->setChildIndex(0, i);
    bottom->setChildReversed(0, false);
    bottom->setCoordinates(i * 10, 10);
    top->setBottomParseIndex(NULL_INDEX);
    top->setParentIndex(i);
    top->setCoordinates(i * 10, 10);
    top->setParentReversed(false);
    top->setNextParalogyIndex(NULL_INDEX);
    bottom->toRight();
    top->toRight();
  }

  // single inversion at segment [50,59]
  bottom = ancGenome->getBottomSegmentIterator(5);
  top->toChild(bottom, 0);
  bottom->setChildReversed(0, true);
  top->setParentReversed(true);
}

void InversionTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);

  const Genome* ancGenome = alignment->openGenome("AncGenome");

  SGBuilder sgBuild;
  
  sgBuild.init(alignment, ancGenome, false, false);
  sgBuild.addGenome(ancGenome);
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  sgBuild.addGenome(leaf1Genome);

  const SideGraph* sg = sgBuild.getSideGraph();

  // only one sequnce in our graph
  CuAssertTrue(_testCase, sg->getNumSequences() == 1);
  const SGSequence* seq = sg->getSequence(0);
  CuAssertTrue(_testCase, seq->getLength() == ancGenome->getSequenceLength());
  CuAssertTrue(_testCase, seq->getName() ==
               ancGenome->getSequenceBySite(0)->getFullName());

  const SideGraph::JoinSet* joins = sg->getJoinSet();
  
  // single inversion at segment [50,59]
  // we expect a join connecting 49-forward to 59-forward
  SGJoin j = SGJoin(SGSide(SGPosition(0, 49), true),
                    SGSide(SGPosition(0, 59), true));
  const SGJoin* join1 = sg->getJoin(&j);
  CuAssertTrue(_testCase, join1 != NULL);
  // and we expect a join connecting 50-rev to 60-rev
  j = SGJoin(SGSide(SGPosition(0, 50), false),
             SGSide(SGPosition(0, 60), false));
  const SGJoin* join2 = sg->getJoin(&j);
  CuAssertTrue(_testCase, join2 != NULL);
  // and nothing else
  CuAssertTrue(_testCase, joins->size() == 2);

  const vector<const Sequence*>& halSequences = sgBuild.getHalSequences();
  for (size_t i = 0; i < halSequences.size(); ++i)
  {
    vector<SGSegment> path;
    sgBuild.getHalSequencePath(halSequences[i], path);
    CuAssertTrue(_testCase, sgBuild.verifyPath(halSequences[i], path) == true);
  }
}

void sgBuilderInversionTest(CuTest *testCase)
{
  try
  {
    InversionTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}


///////////////////////////////////////////////////////////////////////////
//
//           HARDER INVERISON TEST (more inversions more sequences)
//
///////////////////////////////////////////////////////////////////////////

struct Inversion2Test : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

void Inversion2Test::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  Genome* leaf1Genome = alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);

  vector<Sequence::Info> seqVec(2);
  seqVec[0] =Sequence::Info("DummySequence", 30, 0, 3);
  seqVec[1] =Sequence::Info("AncSequence", 100, 0, 10);
  ancGenome->setDimensions(seqVec);

  seqVec[0] =Sequence::Info("HeadSequence", 10, 1, 0);
  seqVec[1] =Sequence::Info("LeafSequence", 100, 10, 0);
  seqVec.push_back(Sequence::Info("TailSequence", 20, 2, 0));
  leaf1Genome->setDimensions(seqVec);

  string dna = randDNA(ancGenome->getSequenceLength());
  ancGenome->setString(dna);
  string leafdna = dna.substr(0, 10);
  string temp = dna.substr(30, 10);
  reverseComplement(temp);
  leafdna += temp;
  leafdna += dna.substr(40, 40);
  temp = dna.substr(80, 10);
  reverseComplement(temp);
  leafdna += temp;
  leafdna += dna.substr(90, 30);
  leafdna += dna.substr(80, 10);
  leafdna += dna.substr(10, 20);
  leaf1Genome->setString(leafdna);

  TopSegmentIteratorPtr top = leaf1Genome->getTopSegmentIterator();
  BottomSegmentIteratorPtr bottom = ancGenome->getBottomSegmentIterator();

  vector<size_t> downMap(ancGenome->getNumBottomSegments());
  vector<size_t> upMap(leaf1Genome->getNumTopSegments());
  assert(downMap.size() == upMap.size());
  downMap[0] = 0;
  upMap[0] = 0;
  downMap[1] = 11;
  upMap[11] = 1;
  downMap[2] = 12;
  upMap[12] = 2;
  for (size_t i = 0; i < 10; ++i)
  {
    downMap[3 + i] = 1 + i;
    upMap[1 + i] = 3 + i;
  }

  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    bottom->setTopParseIndex(NULL_INDEX);
    bottom->setChildIndex(0, downMap[i]);
    bottom->setChildReversed(0, false);
    bottom->setCoordinates(i * 10, 10);
    top->setBottomParseIndex(NULL_INDEX);
    top->setParentIndex(upMap[i]);
    top->setCoordinates(i * 10, 10);
    top->setParentReversed(false);
    top->setNextParalogyIndex(NULL_INDEX);
    bottom->toRight();
    top->toRight();
  }

  // single inversion at segment [50,59] (in ancRefSequence)
  bottom = ancGenome->getBottomSegmentIterator(5 + 3);
  top->toChild(bottom, 0);
  bottom->setChildReversed(0, true);
  top->setParentReversed(true);

  // throw in a dupe / deletion
  hal_index_t topIdx = top->getArrayIndex();
  top->setNextParalogyIndex(top->getArrayIndex() + 4);
  top = leaf1Genome->getTopSegmentIterator(top->getArrayIndex() + 4);
  top->setNextParalogyIndex(topIdx);
  bottom->toParent(top);
  bottom->setChildIndex(0, NULL_INDEX);
  bottom->setChildReversed(0, false);
  top->setParentIndex(5 + 3);
  top->setParentReversed(false);
  
  // inversion at first segment
  bottom = ancGenome->getBottomSegmentIterator(0 + 3);
  top->toChild(bottom, 0);
  bottom->setChildReversed(0, true);
  top->setParentReversed(true);
}

void Inversion2Test::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* ancGenome = alignment->openGenome("AncGenome");

  SGBuilder sgBuild;
  
  sgBuild.init(alignment, ancGenome, false, false);
  sgBuild.addGenome(ancGenome);
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  sgBuild.addGenome(leaf1Genome);

  const SideGraph* sg = sgBuild.getSideGraph();

  // only two sequnce in our graph
  CuAssertTrue(_testCase, sg->getNumSequences() == 2);
  const Sequence* halSeq = ancGenome->getSequenceBySite(50);
  CuAssertTrue(_testCase, halSeq->getName() == "AncSequence");
  const SGSequence* seq = sg->getSequence(1);
  CuAssertTrue(_testCase, seq->getName() == halSeq->getFullName());
  CuAssertTrue(_testCase, seq->getLength() == halSeq->getSequenceLength());

  const SideGraph::JoinSet* joins = sg->getJoinSet();

  size_t joinCount = 0;
  // single inversion at segment [50,59]
  // we expect a join connecting 49-forward to 59-forward
  SGJoin j = SGJoin(SGSide(SGPosition(1, 49), true),
                    SGSide(SGPosition(1, 59), true));
  const SGJoin* join1 = sg->getJoin(&j);
  CuAssertTrue(_testCase, join1 != NULL);
  ++joinCount;
  // and we expect a join connecting 50-rev to 60-rev
  j = SGJoin(SGSide(SGPosition(1, 50), false),
             SGSide(SGPosition(1, 60), false));
  const SGJoin* join2 = sg->getJoin(&j);
  CuAssertTrue(_testCase, join2 != NULL);
  ++joinCount;
  
  // inversion at segment[0, 9]
  j = SGJoin(SGSide(SGPosition(1, 0), false),
             SGSide(SGPosition(1, 10), false));
  join1 = sg->getJoin(&j);
  CuAssertTrue(_testCase, join1 != NULL);
  ++joinCount;
  
  // reverse dupe of first inverion
  j = SGJoin(SGSide(SGPosition(1, 89), true),
             SGSide(SGPosition(1, 50), false));
  const SGJoin* join = sg->getJoin(&j);
  CuAssertTrue(_testCase, join != NULL);
  ++joinCount;
  
  // and nothing else
  CuAssertTrue(_testCase, joins->size() == joinCount);

  const vector<const Sequence*>& halSequences = sgBuild.getHalSequences();
  for (size_t i = 0; i < halSequences.size(); ++i)
  {
    vector<SGSegment> path;
    sgBuild.getHalSequencePath(halSequences[i], path);
    CuAssertTrue(_testCase, sgBuild.verifyPath(halSequences[i], path) == true);
  }
}

void sgBuilderInversion2Test(CuTest *testCase)
{
  try
  {
    Inversion2Test tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}


///////////////////////////////////////////////////////////////////////////
//
//           INSERTION TEST (as added to inversionTest2)
//
///////////////////////////////////////////////////////////////////////////

struct InsertionTest : public Inversion2Test
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

void InsertionTest::createCallBack(AlignmentPtr alignment)
{
  Inversion2Test::createCallBack(alignment);
    
  Genome* ancGenome = alignment->openGenome("AncGenome");
  Genome* leaf1Genome = alignment->openGenome("Leaf1");

  TopSegmentIteratorPtr top = leaf1Genome->getTopSegmentIterator();
  BottomSegmentIteratorPtr bottom = ancGenome->getBottomSegmentIterator();

  //spike in insertion at beginning
  top = leaf1Genome->getTopSegmentIterator(0);
  bottom->toParent(top);
  top->setParentIndex(NULL_INDEX);
  bottom->setChildIndex(0, NULL_INDEX);

  //spike in insertion at beginning of second sequence
  top = leaf1Genome->getTopSegmentIterator(1);
  bottom->toParent(top);
  top->setParentIndex(NULL_INDEX);
  bottom->setChildIndex(0, NULL_INDEX);

  //spike in insertion at middle of second sequence
  top = leaf1Genome->getTopSegmentIterator(3);
  bottom->toParent(top);
  top->setParentIndex(NULL_INDEX);
  bottom->setChildIndex(0, NULL_INDEX);

  //spike in insertion at end
  top = leaf1Genome->getTopSegmentIterator(9);
  bottom->toParent(top);
  top->setParentIndex(NULL_INDEX);
  bottom->setChildIndex(0, NULL_INDEX);

}

void InsertionTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* ancGenome = alignment->openGenome("AncGenome");

  SGBuilder sgBuild;
  
  sgBuild.init(alignment, ancGenome, false, false);
  sgBuild.addGenome(ancGenome);
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  sgBuild.addGenome(leaf1Genome);

  const SideGraph* sg = sgBuild.getSideGraph();

  // only two sequnce in our graph + three inserted. 
  CuAssertTrue(_testCase, sg->getNumSequences() == 6);
  
  const vector<const Sequence*>& halSequences = sgBuild.getHalSequences();
  for (size_t i = 0; i < halSequences.size(); ++i)
  {
    vector<SGSegment> path;
    sgBuild.getHalSequencePath(halSequences[i], path);
    CuAssertTrue(_testCase, sgBuild.verifyPath(halSequences[i], path) == true);
  }

}

void sgBuilderInsertionTest(CuTest *testCase)
{
  try
  {
    InsertionTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

///////////////////////////////////////////////////////////////////////////
//
//           Empty sequence tests
//
///////////////////////////////////////////////////////////////////////////

struct EmptySequenceTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

void EmptySequenceTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  Genome* leaf1Genome = alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);

  vector<Sequence::Info> seqVec(2);
  seqVec[0] =Sequence::Info("DummySequence", 30, 0, 3);
  seqVec[1] =Sequence::Info("AncSequence", 100, 0, 10);
  seqVec.push_back(Sequence::Info("EmtpyRootSequence", 0, 0, 0));
  ancGenome->setDimensions(seqVec);

  seqVec[0] =Sequence::Info("HeadSequence", 10, 1, 0);
  seqVec[1] =Sequence::Info("LeafSequence", 100, 10, 0);
  seqVec.push_back(Sequence::Info("EmtpyLeafSequence", 0, 0, 0));
  seqVec.push_back(Sequence::Info("TailSequence", 20, 2, 0));
  leaf1Genome->setDimensions(seqVec);

  string dna = randDNA(ancGenome->getSequenceLength());
  ancGenome->setString(dna);
  leaf1Genome->setString(dna);

  TopSegmentIteratorPtr top = leaf1Genome->getTopSegmentIterator();
  BottomSegmentIteratorPtr bottom = ancGenome->getBottomSegmentIterator();

  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    bottom->setTopParseIndex(NULL_INDEX);
    bottom->setChildIndex(0, i);
    bottom->setChildReversed(0, false);
    bottom->setCoordinates(i * 10, 10);
    top->setBottomParseIndex(NULL_INDEX);
    top->setParentIndex(i);
    top->setCoordinates(i * 10, 10);
    top->setParentReversed(false);
    top->setNextParalogyIndex(NULL_INDEX);
    bottom->toRight();
    top->toRight();
  }
}


void EmptySequenceTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* ancGenome = alignment->openGenome("AncGenome");

  SGBuilder sgBuild;
  
  sgBuild.init(alignment, ancGenome, false, false);
  sgBuild.addGenome(ancGenome);
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  sgBuild.addGenome(leaf1Genome);

  const SideGraph* sg = sgBuild.getSideGraph();

  // we expect empty sequences to not be converted.  
  CuAssertTrue(_testCase, sg->getNumSequences() == 2);

  const SideGraph::JoinSet* joins = sg->getJoinSet();

  // no breakpoints except sequence size difference
  CuAssertTrue(_testCase, joins->size() == 1);

  const vector<const Sequence*>& halSequences = sgBuild.getHalSequences();
  for (size_t i = 0; i < halSequences.size(); ++i)
  {
    vector<SGSegment> path;
    sgBuild.getHalSequencePath(halSequences[i], path);
    CuAssertTrue(_testCase, sgBuild.verifyPath(halSequences[i], path) == true);
  }

}

void sgBuilderEmptySequenceTest(CuTest *testCase)
{
  try
  {
    EmptySequenceTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

///////////////////////////////////////////////////////////////////////////
//
//            SNP TEST
//
///////////////////////////////////////////////////////////////////////////

struct SNPTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

void SNPTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  Genome* leaf1Genome = alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);

  vector<Sequence::Info> seqVec(1);
  seqVec[0] =Sequence::Info("AncSequence", 100, 0, 10);
  ancGenome->setDimensions(seqVec);

  seqVec[0] =Sequence::Info("LeafSequence", 100, 10, 0);
  leaf1Genome->setDimensions(seqVec);

  string dna = randDNA(ancGenome->getSequenceLength());
  string leafdna = dna;
  
  // Spike in point mutation from A->C at position 5
  dna[5] = 'A';
  leafdna[5] = 'C';
  
  ancGenome->setString(dna);
  leaf1Genome->setString(leafdna);

  TopSegmentIteratorPtr top = leaf1Genome->getTopSegmentIterator();
  BottomSegmentIteratorPtr bottom = ancGenome->getBottomSegmentIterator();

  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    bottom->setTopParseIndex(NULL_INDEX);
    bottom->setChildIndex(0, i);
    bottom->setChildReversed(0, false);
    bottom->setCoordinates(i * 10, 10);
    top->setBottomParseIndex(NULL_INDEX);
    top->setParentIndex(i);
    top->setCoordinates(i * 10, 10);
    top->setParentReversed(false);
    top->setNextParalogyIndex(NULL_INDEX);
    bottom->toRight();
    top->toRight();
  }
}

void SNPTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);

  const Genome* ancGenome = alignment->openGenome("AncGenome");

  SGBuilder sgBuild;
  
  sgBuild.init(alignment, ancGenome, false, false);
  sgBuild.addGenome(ancGenome);
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  sgBuild.addGenome(leaf1Genome);

  const SideGraph* sg = sgBuild.getSideGraph();

  const vector<const Sequence*>& halSequences = sgBuild.getHalSequences();
  for (size_t i = 0; i < halSequences.size(); ++i)
  {
    vector<SGSegment> path;
    sgBuild.getHalSequencePath(halSequences[i], path);
    CuAssertTrue(_testCase, sgBuild.verifyPath(halSequences[i], path) == true);
  }

  CuAssertTrue(_testCase, sg->getNumSequences() == 2);

  SGJoin join;
  join = SGJoin(SGSide(SGPosition(0, 4), true),
                SGSide(SGPosition(1, 0), false));
  CuAssertTrue(_testCase, sg->getJoin(&join) != NULL);

  join = SGJoin(SGSide(SGPosition(1, 0), true),
                SGSide(SGPosition(0, 6), false));
  CuAssertTrue(_testCase, sg->getJoin(&join) != NULL);

  CuAssertTrue(_testCase, sg->getJoinSet()->size() == 2);
}

void sgBuilderSNPTest(CuTest *testCase)
{
  try
  {
    SNPTest tester;
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
  SUITE_ADD_TEST(suite, sgBuilderInversionTest);
  SUITE_ADD_TEST(suite, sgBuilderInversion2Test);
  SUITE_ADD_TEST(suite, sgBuilderInsertionTest);
  SUITE_ADD_TEST(suite, sgBuilderEmptySequenceTest);
  SUITE_ADD_TEST(suite, sgBuilderSNPTest);
  return suite;
}
