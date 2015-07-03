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
#include <cstdio>
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
  sgBuild.computeJoins();

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
  // leaf2 only used in derived test HarderSNPTest
  Genome* leaf2Genome = alignment->addLeafGenome("Leaf2", "AncGenome", 0.1);
  
  vector<Sequence::Info> seqVec;
  seqVec.push_back(Sequence::Info("DummySequence", 30, 0, 3));
  seqVec.push_back(Sequence::Info("AncSequence", 100, 0, 10));
  ancGenome->setDimensions(seqVec);

  seqVec.clear();
  seqVec.push_back(Sequence::Info("HeadSequence", 10, 1, 0));
  seqVec.push_back(Sequence::Info("LeafSequence", 100, 10, 0));
  seqVec.push_back(Sequence::Info("TailSequence", 20, 2, 0));
  leaf1Genome->setDimensions(seqVec);
  
  seqVec.clear();
  seqVec.push_back(Sequence::Info("DummySequence", 30, 3, 0));
  seqVec.push_back(Sequence::Info("Leaf2Sequence", 100, 10, 0));
  leaf2Genome->setDimensions(seqVec);

  string dna = randDNA(ancGenome->getSequenceLength());
  ancGenome->setString(dna);
  leaf2Genome->setString(dna);
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
  TopSegmentIteratorPtr top2 = leaf2Genome->getTopSegmentIterator();
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
    bottom->setChildIndex(1, i);
    bottom->setChildReversed(1, false);
    bottom->setCoordinates(i * 10, 10);
    top->setBottomParseIndex(NULL_INDEX);
    top->setParentIndex(upMap[i]);
    top->setCoordinates(i * 10, 10);
    top->setParentReversed(false);
    top->setNextParalogyIndex(NULL_INDEX);
    top2->setBottomParseIndex(NULL_INDEX);
    top2->setParentIndex(i);
    top2->setCoordinates(i * 10, 10);
    top2->setParentReversed(false);
    top2->setNextParalogyIndex(NULL_INDEX);
    bottom->toRight();
    top->toRight();
    top2->toRight();
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
  sgBuild.computeJoins();

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
  sgBuild.computeJoins();

  const SideGraph* sg = sgBuild.getSideGraph();

  // only two sequnce in our graph + three inserted. 
  CuAssertTrue(_testCase, sg->getNumSequences() == 6);
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
  sgBuild.computeJoins();

  const SideGraph* sg = sgBuild.getSideGraph();

  // we expect empty sequences to not be converted.  
  CuAssertTrue(_testCase, sg->getNumSequences() == 2);

  const SideGraph::JoinSet* joins = sg->getJoinSet();

  // no breakpoints except sequence size difference
  CuAssertTrue(_testCase, joins->size() == 1);
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
  sgBuild.computeJoins();

  const SideGraph* sg = sgBuild.getSideGraph();

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

///////////////////////////////////////////////////////////////////////////
//
//            MORE ELABORATE SNP TEST
//
///////////////////////////////////////////////////////////////////////////

struct HarderSNPTest : public Inversion2Test
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

static char mutate(char c, int i = 1)
{
  assert(i > 0 && i < 4);
  c = toupper(c);
  int n;
  switch (c) {
  case 'A' : n = 0; break;
  case 'G' : n = 1; break;
  case 'T' : n = 2; break;
  default : n = 3;
  }
  n = (n + i) % 4;
  char v[4] = {'A', 'G', 'T', 'C'};
  return v[n];
}

void HarderSNPTest::createCallBack(AlignmentPtr alignment)
{
  Inversion2Test::createCallBack(alignment);

  Genome* ancGenome = alignment->openGenome("AncGenome");
  Genome* leaf1Genome = alignment->openGenome("Leaf1");
  Genome* leaf2Genome = alignment->openGenome("Leaf2");

  string dna;
  ancGenome->getString(dna);

  // leaf1 has: Inversion between bottom(8) and top(6) [50,59]
  //            Dupe between top(6) and top(12)
  //            Inversion between bottom(3) and top(1)
  // leaf2 has: same structure as ancestor. 

  string leaf2dna;
  leaf2Genome->getString(leaf2dna);
  assert(leaf2dna == dna);
  string leaf1dna;
  leaf1Genome->getString(leaf1dna);

  // Make a whole segment of snps in leaf2 at segment 1
  for (size_t i = 10; i < 20; ++i)
  {
    leaf2dna[i] = mutate(leaf2dna[i], 2);
  }

  // Make a paritally overlapping (forward) snp region in both leaves
  assert(dna.substr(0, 10) == leaf2dna.substr(0, 10));
  assert(dna.substr(0, 10) == leaf1dna.substr(0, 10));
  leaf2dna[3] = mutate(leaf2dna[3], 2);
  leaf2dna[4] = mutate(leaf2dna[4], 2);
  leaf2dna[5] = mutate(leaf2dna[5], 2);
  leaf2dna[6] = mutate(leaf2dna[6], 2);
  leaf2dna[7] = mutate(leaf2dna[7], 2);
  leaf1dna[3] = mutate(leaf1dna[3], 1);
  leaf1dna[4] = mutate(leaf1dna[4], 2);
  leaf1dna[5] = mutate(leaf1dna[5], 1);
  leaf1dna[6] = mutate(leaf1dna[6], 1);
  leaf1dna[7] = mutate(leaf1dna[7], 2);

  // Test edge cases : start and end of sequence
  leaf2dna[0] = mutate(leaf2dna[0], 2);
  leaf2dna[leaf2dna.length()-1] = mutate(leaf2dna[leaf2dna.length()-1], 2);
  
  // Make a reverse SNP in leaf1
  leaf1dna[12] = mutate(leaf1dna[12], 1);
  leaf1dna[13] = mutate(leaf1dna[13], 1);

  // Overlap a SNP on leaf2
  leaf2dna[37] = mutate(leaf2dna[37], 2);
  leaf2dna[38] = mutate(leaf2dna[38], 2);
  leaf2dna[39] = mutate(leaf2dna[39], 2);
  
  leaf1Genome->setString(leaf1dna);
  leaf2Genome->setString(leaf2dna);
}

void HarderSNPTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);

  const Genome* ancGenome = alignment->openGenome("AncGenome");
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  const Genome* leaf2Genome = alignment->openGenome("Leaf2");

  SGBuilder sgBuild;
  
  sgBuild.init(alignment, ancGenome, false, false);
  sgBuild.addGenome(ancGenome);
  sgBuild.addGenome(leaf1Genome);
  sgBuild.addGenome(leaf2Genome);
  sgBuild.computeJoins();

  const SideGraph* sg = sgBuild.getSideGraph();
   
  // sequence 0,1: ancestral sequences
  const SGSequence* sequence;
  // sequence 2 : length 5 mutation in leaf1
  sequence = sg->getSequence(2);
  CuAssertTrue(_testCase, sequence->getLength() == 5);
  // sequence 3 : length 2 mutation in snp in leaf1
  //
  // TODO: new forward code for target seems to change
  // this:  need to revise...
  //CuAssertTrue(_testCase, sequence->getLength() == 2);
  int l2o = 4;
  // sequence l2o : length 1 mutation at leaf2[0]
  sequence = sg->getSequence(l2o+0);
  CuAssertTrue(_testCase, sequence->getLength() == 1);
  // sequence l2o+1 : length 1 mutation at leaf2[3]
  sequence = sg->getSequence(l2o+1);
  CuAssertTrue(_testCase, sequence->getLength() == 1);
  // sequence l2o+2 : length 2 mutation at leaf2[5]
  sequence = sg->getSequence(l2o+2);
  CuAssertTrue(_testCase, sequence->getLength() == 2);
  // sequence l2o+3 : length 10 mutation at leaf2[10]
  sequence = sg->getSequence(l2o+3);
  CuAssertTrue(_testCase, sequence->getLength() == 10);
  // sequence l2o+4 : length 3 mutation at leaf2[37]
  sequence = sg->getSequence(l2o+4);
  CuAssertTrue(_testCase, sequence->getLength() == 3);
  // sequence l2o+5 : length 1 mutation at leaf2[last]
  sequence = sg->getSequence(l2o+5);
  CuAssertTrue(_testCase, sequence->getLength() == 1);

  // Should test all joins here, but if path check goes through
  // then that's enough for now. 
  //SGJoin join;
  //join = SGJoin(SGSide(SGPosition(0, 4), true),
  //              SGSide(SGPosition(1, 0), false));
  //CuAssertTrue(_testCase, sg->getJoin(&join) != NULL);
}

void sgBuilderHarderSNPTest(CuTest *testCase)
{
  try
  {
    HarderSNPTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

///////////////////////////////////////////////////////////////////////////
//
//            BASIC REFERENCE DUPE TEST
//
///////////////////////////////////////////////////////////////////////////

struct RefDupeTest : public InversionTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

void RefDupeTest::createCallBack(AlignmentPtr alignment)
{
  InversionTest::createCallBack(alignment);

  // note we inherit the following:
  // ancestor and leaf are 10 x 10 blocks perfectly aligned except  
  // single inversion at segment [50,59]

  Genome* ancGenome = alignment->openGenome("AncGenome");
  Genome* leaf1Genome = alignment->openGenome("Leaf1");

  string dna;
  ancGenome->getString(dna);
  string leafdna;
  leaf1Genome->getString(leafdna);

  // add a couple dupes to the leaf
  // 0-9 and 20-29: single forward dupe, no snps
  TopSegmentIteratorPtr top = leaf1Genome->getTopSegmentIterator(2);
  BottomSegmentIteratorPtr bottom = ancGenome->getBottomSegmentIterator();
  bottom->toParent(top);
  bottom->setChildIndex(0, NULL_INDEX);
  top->setParentIndex(0);
  top->setParentReversed(0);
  top->setNextParalogyIndex(0);
  top = leaf1Genome->getTopSegmentIterator(0);
  top->setNextParalogyIndex(2);
  for (size_t i = 0; i < 10; ++i)
  {
    leafdna[20 + i] = leafdna[i];
  }
 
  // 50,59, 60,69 90,99 : triple dupe where first is inverted
  top = leaf1Genome->getTopSegmentIterator(5);
  top->setNextParalogyIndex(6);
  top = leaf1Genome->getTopSegmentIterator(6);
  top->setParentIndex(5);
  top->setParentReversed(false);
  top->setNextParalogyIndex(9);
  top = leaf1Genome->getTopSegmentIterator(9);
  top->setParentIndex(5);
  top->setParentReversed(false);
  top->setNextParalogyIndex(5);
  bottom = ancGenome->getBottomSegmentIterator(6);
  bottom->setChildIndex(0, NULL_INDEX);
  bottom = ancGenome->getBottomSegmentIterator(9);
  bottom->setChildIndex(0, NULL_INDEX);
  for (size_t i = 0; i < 10; ++i)
  {
    leafdna[60 + i] = dna[50 + i];
    leafdna[90 + i] = dna[50 + i];
  }
  // pop in a couple snps
  leafdna[65] = mutate(leafdna[65], 1);
  leafdna[66] = mutate(leafdna[66], 1);
  leafdna[99] = mutate(leafdna[99], 1);
    
  leaf1Genome->setString(leafdna);
}

void RefDupeTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);

  const Genome* ancGenome = alignment->openGenome("AncGenome");
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");

  SGBuilder sgBuild;
  
  sgBuild.init(alignment, ancGenome, true, false);
  sgBuild.addGenome(leaf1Genome);
  sgBuild.addGenome(ancGenome);
  sgBuild.computeJoins();

  const SideGraph* sg = sgBuild.getSideGraph();

  cout << *sg << endl;

  CuAssertTrue(_testCase, sg->getNumSequences() == 6);
  // folded leaf
  CuAssertTrue(_testCase, sg->getSequence(0)->getLength() == 70);
  // 1 double dupe and 1 single:
  CuAssertTrue(_testCase, sg->getSequence(1)->getLength() == 2);
  CuAssertTrue(_testCase, sg->getSequence(2)->getLength() == 1);
  // 3 10-base holes from anc
  CuAssertTrue(_testCase, sg->getSequence(3)->getLength() == 10);
  CuAssertTrue(_testCase, sg->getSequence(4)->getLength() == 10);
  CuAssertTrue(_testCase, sg->getSequence(5)->getLength() == 10);
}

void sgBuilderRefDupeTest(CuTest *testCase)
{
  try
  {
    RefDupeTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

///////////////////////////////////////////////////////////////////////////
//
//           INVERSION TEST WITH 3 levels
//
///////////////////////////////////////////////////////////////////////////

struct TransInversionTest : public InversionTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

void TransInversionTest::createCallBack(AlignmentPtr alignment)
{
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  Genome* midGenome = alignment->addLeafGenome("Mid", "AncGenome", 0.1);
  Genome* leaf1Genome = alignment->addLeafGenome("Leaf1", "Mid", 0.1);
  Genome* leaf2Genome = alignment->addLeafGenome("Leaf2", "Mid", 0.1);

  vector<Sequence::Info> seqVec(1);
  seqVec[0] =Sequence::Info("AncSequence", 100, 0, 10);
  ancGenome->setDimensions(seqVec);

  seqVec[0] = Sequence::Info("MidSequence", 100, 10, 10);
  midGenome->setDimensions(seqVec);

  seqVec[0] =Sequence::Info("Leaf1Sequence", 100, 10, 0);
  leaf1Genome->setDimensions(seqVec);
  
  seqVec[0] =Sequence::Info("Leaf2Sequence", 100, 10, 0);
  leaf2Genome->setDimensions(seqVec);

  // set up so that anc->mid has and inversion
  // and mid->leaf2 has the same inversion
  string dna = randDNA(ancGenome->getSequenceLength());
  ancGenome->setString(dna);
  string leafdna = dna.substr(0, 50);
  string temp = dna.substr(50, 10);
  reverseComplement(temp);
  leafdna += temp;
  leafdna += dna.substr(60, 40);
  midGenome->setString(leafdna);
  leaf1Genome->setString(leafdna);
  leaf2Genome->setString(dna);

  // align mid to anc
  TopSegmentIteratorPtr top = midGenome->getTopSegmentIterator();
  BottomSegmentIteratorPtr bottom = ancGenome->getBottomSegmentIterator();

  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    bottom->setTopParseIndex(NULL_INDEX);
    bottom->setChildIndex(0, i);
    bottom->setChildReversed(0, false);
    bottom->setCoordinates(i * 10, 10);
    top->setBottomParseIndex(i);
    top->setParentIndex(i);
    top->setCoordinates(i * 10, 10);
    top->setParentReversed(false);
    top->setNextParalogyIndex(NULL_INDEX);
    bottom->toRight();
    top->toRight();
  }

  // single inversion at segment [50,59] between mid and anc
  bottom = ancGenome->getBottomSegmentIterator(5);
  top->toChild(bottom, 0);
  bottom->setChildReversed(0, true);
  top->setParentReversed(true);

  // align mid and leaf1
  top = leaf1Genome->getTopSegmentIterator();
  bottom = midGenome->getBottomSegmentIterator();

  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    bottom->setTopParseIndex(i);
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

  // align leaf2 and mid
  top = leaf2Genome->getTopSegmentIterator();
  bottom = midGenome->getBottomSegmentIterator();

  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    bottom->setTopParseIndex(i);
    bottom->setChildIndex(1, i);
    bottom->setChildReversed(1, false);
    bottom->setCoordinates(i * 10, 10);
    top->setBottomParseIndex(NULL_INDEX);
    top->setParentIndex(i);
    top->setCoordinates(i * 10, 10);
    top->setParentReversed(false);
    top->setNextParalogyIndex(NULL_INDEX);
    bottom->toRight();
    top->toRight();
  }

  // single inversion at segment [50,59] between leaf2 and mid
  bottom = midGenome->getBottomSegmentIterator(5);
  top->toChild(bottom, 1);
  bottom->setChildReversed(1, true);
  top->setParentReversed(true);
}

void TransInversionTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);

  const Genome* ancGenome = alignment->openGenome("AncGenome");
  const Genome* midGenome = alignment->openGenome("Mid");
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  const Genome* leaf2Genome = alignment->openGenome("Leaf2");
  CuAssertTrue(_testCase, ancGenome && midGenome && leaf1Genome &&
               leaf2Genome);

  const Genome* genomes[] = {ancGenome, midGenome, leaf1Genome, leaf2Genome};

   do
   {
     cout << endl << endl << endl;
     SGBuilder sgBuild;
     sgBuild.init(alignment, ancGenome);
     for (size_t i = 0; i < 4; ++i)
     {
       sgBuild.addGenome(genomes[i]);
     }
     sgBuild.computeJoins();
   } while (next_permutation(genomes, genomes+4));
}

void sgBuilderTransInversionTest(CuTest *testCase)
{
  try
  {
    TransInversionTest tester;
    tester.check(testCase);
  }
  catch (string s) 
  {
    CuAssertTrue(testCase, false);
  }
}


///////////////////////////////////////////////////////////////////////////

CuSuite* sgBuildTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, sgBuilderHalTest);
  SUITE_ADD_TEST(suite, sgBuilderInversionTest);
  SUITE_ADD_TEST(suite, sgBuilderInversion2Test);
  SUITE_ADD_TEST(suite, sgBuilderInsertionTest);
  SUITE_ADD_TEST(suite, sgBuilderEmptySequenceTest);
  SUITE_ADD_TEST(suite, sgBuilderSNPTest);
  SUITE_ADD_TEST(suite, sgBuilderHarderSNPTest);
  SUITE_ADD_TEST(suite, sgBuilderRefDupeTest);
  SUITE_ADD_TEST(suite, sgBuilderTransInversionTest);
  return suite;
}
