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

// Basic test of underlying map structure.  
void snpMapTest(CuTest *testCase)
{
  SGPosition p1(0, 10);
  SGPosition p2(0, 1);
  SGPosition p3(1, 10);
  SGPosition p4(1, 0);
  SGPosition p5(1, 3);
  SGPosition p6(2, 0);
  SGPosition p7(2, 1);
  
  SNPHandler snpHandler(NULL, false);

  snpHandler.addSNP(p1, 'A', p4);
  snpHandler.addSNP(p3, 'A', p2);
  snpHandler.addSNP(p1, 'c', p6);
  snpHandler.addSNP(p1, 'T', p5);

  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'A') == p4);
  CuAssertTrue(testCase, snpHandler.findSNP(p3, 'A') == p2);
  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'c') == p6);
  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'C') == p6);
  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'G') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandler.findSNP(p3, 'G') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandler.findSNP(p2, 'A') == p2);

  SNPHandler snpHandlerCS(NULL, true);

  snpHandlerCS.addSNP(p1, 'A', p4);
  snpHandlerCS.addSNP(p3, 'A', p2);
  snpHandlerCS.addSNP(p1, 'c', p6);
  snpHandlerCS.addSNP(p1, 'T', p5);
  snpHandlerCS.addSNP(p1, 't', p7);

  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'A') == p4);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p3, 'A') == p2);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'c') == p6);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'C') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'T') == p5);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 't') == p7);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'G') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p3, 'G') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p2, 'A') == p2);
}

/** easiest case: we add a single SNP
 */
void snpHandlerSingleSNPTest(CuTest *testCase)
{
  SideGraph sg;
  sg.addSequence(new SGSequence(-1, 100, "Seq0"));
  SGLookup lookup;
  vector<string> seqNames;
  seqNames.push_back("Seq0");
  seqNames.push_back("Seq1");  
  lookup.init(seqNames);
  SNPHandler snpHandler(&sg);
  
  string dna = "A";
  SGPosition srcPos(1, 5);
  SGPosition tgtPos(0, 5);

  pair<SGSide, SGSide> hooks = snpHandler.createSNP(dna, 0, 1, srcPos, tgtPos,
                                                    false, &lookup);

  cout << "outhooks " << hooks.first << " , " << hooks.second << endl;
  
  CuAssertTrue(testCase, hooks.first.getBase() == hooks.second.getBase());
  CuAssertTrue(testCase, hooks.first.getForward() == false);
  CuAssertTrue(testCase, hooks.second.getForward() == true);
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos, 'A') ==
               hooks.first.getBase());
  CuAssertTrue(testCase, sg.getNumSequences() == 2);
  CuAssertTrue(testCase, hooks.first.getBase().getSeqID() == 1);
  CuAssertTrue(testCase, sg.getSequence(1)->getLength() == 1);
  CuAssertTrue(testCase, lookup.mapPosition(srcPos) == hooks.second);
}

/** next case: we add a single multibase SNP
 */
void snpHandlerMultibaseSNPTest(CuTest *testCase)
{
  SideGraph sg;
  sg.addSequence(new SGSequence(-1, 100, "Seq0"));
  SGLookup lookup;
  vector<string> seqNames;
  seqNames.push_back("Seq0");
  seqNames.push_back("Seq1");  
  lookup.init(seqNames);
  SNPHandler snpHandler(&sg);
  
  string dna = "ACCCTA";
  SGPosition srcPos(1, 5);
  SGPosition srcPos2(1, 7);
  SGPosition tgtPos(0, 5);
  SGPosition tgtPos1(0, 6);
  SGPosition tgtPos2(0, 7);

  pair<SGSide, SGSide> hooks = snpHandler.createSNP(dna, 2, 3, srcPos, tgtPos,
                                                    false, &lookup);

  cout << "outhooks " << hooks.first << " , " << hooks.second << endl;
  
  CuAssertTrue(testCase, hooks.first.getForward() == false);
  CuAssertTrue(testCase, hooks.second.getForward() == true);
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos, 'C') ==
               hooks.first.getBase());
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos2, 'T') ==
               hooks.second.getBase());
  CuAssertTrue(testCase, sg.getNumSequences() == 2);
  CuAssertTrue(testCase, hooks.first.getBase().getSeqID() == 1);
  CuAssertTrue(testCase, sg.getSequence(1)->getLength() == 3);
  CuAssertTrue(testCase, lookup.mapPosition(srcPos2) == hooks.second);
}

/** next case: we add a single multibase SNP
 */
void snpHandlerMultibaseSNPTest(CuTest *testCase)
{
  SideGraph sg;
  sg.addSequence(new SGSequence(-1, 100, "Seq0"));
  SGLookup lookup;
  vector<string> seqNames;
  seqNames.push_back("Seq0");
  seqNames.push_back("Seq1");  
  lookup.init(seqNames);
  SNPHandler snpHandler(&sg);
  
  string dna = "ACCCTA";
  SGPosition srcPos(1, 5);
  SGPosition srcPos2(1, 7);
  SGPosition tgtPos(0, 5);
  SGPosition tgtPos1(0, 6);
  SGPosition tgtPos2(0, 7);

  pair<SGSide, SGSide> hooks = snpHandler.createSNP(dna, 2, 3, srcPos, tgtPos,
                                                    false, &lookup);

  cout << "outhooks " << hooks.first << " , " << hooks.second << endl;
  
  CuAssertTrue(testCase, hooks.first.getForward() == false);
  CuAssertTrue(testCase, hooks.second.getForward() == true);
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos, 'C') ==
               hooks.first.getBase());
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos2, 'T') ==
               hooks.second.getBase());
  CuAssertTrue(testCase, sg.getNumSequences() == 2);
  CuAssertTrue(testCase, hooks.first.getBase().getSeqID() == 1);
  CuAssertTrue(testCase, sg.getSequence(1)->getLength() == 3);
  CuAssertTrue(testCase, lookup.mapPosition(srcPos2) == hooks.second);
}


CuSuite* snpHandlerTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, snpMapTest);
  SUITE_ADD_TEST(suite, snpHandlerSingleSNPTest);
  SUITE_ADD_TEST(suite, snpHandlerMultibaseSNPTest);
  return suite;
}
