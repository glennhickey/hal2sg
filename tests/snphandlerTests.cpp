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
  const hal::Sequence* halSeq = NULL;
  vector<string> seqNames;
  seqNames.push_back("Seq0");
  seqNames.push_back("Seq1");  
  lookup.init(seqNames);
  SNPHandler snpHandler(&sg);
  
  string srcDNA = "A";
  string tgtDNA = "T";
  SGPosition srcPos(1, 5);
  SGPosition tgtPos(0, 5);

  pair<SGSide, SGSide> hooks = snpHandler.createSNP(srcDNA, tgtDNA, 0, 1,
                                                    halSeq, srcPos, tgtPos,
                                                    false, &lookup, NULL);

  CuAssertTrue(testCase, hooks.first.getBase() == hooks.second.getBase());
  CuAssertTrue(testCase, hooks.first.getForward() == false);
  CuAssertTrue(testCase, hooks.second.getForward() == true);
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos, 'A') ==
               hooks.first.getBase());
  CuAssertTrue(testCase, sg.getNumSequences() == 2);
  CuAssertTrue(testCase, hooks.first.getBase().getSeqID() == 1);
  CuAssertTrue(testCase, sg.getSequence(1)->getLength() == 1);
  CuAssertTrue(testCase, lookup.mapPosition(srcPos) == hooks.second);
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos, 'T') == tgtPos);
}

/** next case: we add a single multibase SNP
 */
void snpHandlerMultibaseSNPTest(CuTest *testCase)
{
  SideGraph sg;
  sg.addSequence(new SGSequence(-1, 100, "Seq0"));
  SGLookup lookup;
  const hal::Sequence* halSeq = NULL;
  vector<string> seqNames;
  seqNames.push_back("Seq0");
  seqNames.push_back("Seq1");  
  lookup.init(seqNames);
  SNPHandler snpHandler(&sg);
  
  string srcDNA = "ACCCTA";
  string tgtDNA = "TGTGGT";
  SGPosition srcPos(1, 5);
  SGPosition srcPos2(1, 7);
  SGPosition tgtPos(0, 5);
  SGPosition tgtPos1(0, 6);
  SGPosition tgtPos2(0, 7);

  pair<SGSide, SGSide> hooks = snpHandler.createSNP(srcDNA, tgtDNA, 2, 3,
                                                    halSeq,
                                                    srcPos, tgtPos,
                                                    false, &lookup, NULL);

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

  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos, 'T') == tgtPos);
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos1, 'G') == tgtPos1);
  CuAssertTrue(testCase, snpHandler.findSNP(tgtPos2, 'G') == tgtPos2);
    
}

/** next case: we add a single multibase SNP and mix in some overlaps
 */
void snpHandlerOverlapSNPTest(CuTest *tc)
{
  SideGraph sg;
  sg.addSequence(new SGSequence(-1, 100, "Seq0"));
  SGLookup lookup;
  const hal::Sequence* halSeq = NULL;
  vector<string> seqNames;
  seqNames.push_back("Seq0");
  seqNames.push_back("Seq1");
  seqNames.push_back("Seq2");
  seqNames.push_back("Seq3");  
  lookup.init(seqNames);
  SNPHandler snpHandler(&sg);

  // this will be our test case here: adding snps from three sequences
  // from the bottom row up:
  string dna0 = "AC...T...CCC.";
  string dna1 = "T...ATA....TT";
  string dna2 = "TC....AGGGGGA";
  string dnaTgt = "GGGGGGGTTTTAG";
  
  SGPosition srcPos(1, 0);
  SGPosition tgtPos(0, 0);
  snpHandler.createSNP(dna0, dnaTgt, 0, 2, halSeq, srcPos, tgtPos, false,
                       &lookup, NULL);
  srcPos.setPos(5);
  tgtPos.setPos(5);
  snpHandler.createSNP(dna0, dnaTgt, 5, 1, halSeq, srcPos, tgtPos, false,
                       &lookup, NULL);
  srcPos.setPos(9);
  tgtPos.setPos(9);
  snpHandler.createSNP(dna0, dnaTgt, 9, 3, halSeq, srcPos, tgtPos, false,
                       &lookup, NULL);

  srcPos = SGPosition(2, 0);
  tgtPos = SGPosition(0, 0);
  snpHandler.createSNP(dna1, dnaTgt, 0, 1, halSeq, srcPos, tgtPos, false,
                       &lookup, NULL);
  srcPos.setPos(4);
  tgtPos.setPos(4);
  snpHandler.createSNP(dna1, dnaTgt, 4, 3, halSeq, srcPos, tgtPos, false,
                       &lookup, NULL);
  srcPos.setPos(11);
  tgtPos.setPos(11);
  snpHandler.createSNP(dna1, dnaTgt, 11, 2, halSeq, srcPos, tgtPos, false,
                       &lookup, NULL);

  srcPos = SGPosition(3, 0);
  tgtPos = SGPosition(0, 0);
  snpHandler.createSNP(dna2, dnaTgt, 0, 2, halSeq, srcPos, tgtPos, false,
                       &lookup, NULL);
  srcPos.setPos(6);
  tgtPos.setPos(6);
  snpHandler.createSNP(dna2, dnaTgt, 6, 7, halSeq, srcPos, tgtPos, false,
                       &lookup, NULL);

  // the snps were added in this order (x's for duplicates).  so we
  // expect the sequence ids that were created to correspond to these
  // numbers
  // "AC...T...CCC.";
  // "T...ATA....TT";
  // "TC....AGGGGGA";
  //
  // 11...2...333.
  // 4...5x6....77
  // xx....x888888

  // top row
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,0), 'A').getSeqID() == 1);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,1), 'C').getSeqID() == 1);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,5), 'T').getSeqID() == 2);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,9), 'C').getSeqID() == 3);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,10), 'C').getSeqID() == 3);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,11), 'C').getSeqID() == 3);

  //middle row
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,0), 'T').getSeqID() == 4);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,4), 'A').getSeqID() == 5);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,5), 'T').getSeqID() == 2);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,6), 'A').getSeqID() == 6);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,11), 'T').getSeqID() == 7);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,12), 'T').getSeqID() == 7);

  //bottom row
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,0), 'T').getSeqID() == 4);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,1), 'C').getSeqID() == 1);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,6), 'A').getSeqID() == 6);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,7), 'G').getSeqID() == 8);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,8), 'G').getSeqID() == 8);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,9), 'G').getSeqID() == 8);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,10), 'G').getSeqID() == 8);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,11), 'G').getSeqID() == 8);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,12), 'A').getSeqID() == 8);

  // spot check position:
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,12), 'A').getPos() == 5);
}

/** next case: we test reverse complement mapping
 */
void snpHandlerInversionSNPTest(CuTest *tc)
{
  SideGraph sg;
  sg.addSequence(new SGSequence(-1, 100, "Seq0"));
  SGLookup lookup;
  const hal::Sequence* halSeq = NULL;
  vector<string> seqNames;
  seqNames.push_back("Seq0");
  seqNames.push_back("Seq1");
  seqNames.push_back("Seq2");
  seqNames.push_back("Seq3");  
  lookup.init(seqNames);
  SNPHandler snpHandler(&sg);

  // single SNP reverse map
  SGPosition srcPos(1, 50);
  SGPosition tgtPos(0, 50);
  std::pair<SGSide, SGSide> ret1 = snpHandler.createSNP("A", "G", 0, 1, halSeq,
                                                        srcPos, tgtPos, true,
                                                        &lookup, NULL);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,50), 'A').getSeqID() == 1);
  // we added the reverse comp. of G because snps for target always forward
  // for now
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,50), 'C').getSeqID() == 0);
  // note: the snp sequence added is walked in the forward direction
  // maybe something to change. 
  CuAssertTrue(tc, ret1.first == SGSide(SGPosition(1, 0), false));
  CuAssertTrue(tc, ret1.second == SGSide(SGPosition(1, 0), true));

  // multi SNP reverse map
  srcPos.setPos(60);
  tgtPos.setPos(60);
  std::pair<SGSide, SGSide> ret2 = snpHandler.createSNP("AACTTC", "AAACGG", 2,
                                                        3, halSeq,
                                                        srcPos, tgtPos, true,
                                                        &lookup, NULL);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,60), 'C').getSeqID() == 2);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,59), 'T').getSeqID() == 2);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,58), 'T').getSeqID() == 2);

  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,60), 'T').getSeqID() == 0);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,59), 'G').getSeqID() == 0);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,58), 'C').getSeqID() == 0);  
  CuAssertTrue(tc, ret2.first == SGSide(SGPosition(2, 0), false));
  CuAssertTrue(tc, ret2.second == SGSide(SGPosition(2, 2), true));

  
  // todo: add overlap cases
  srcPos.setPos(62);
  tgtPos.setPos(58);
//  std::pair<SGSide, SGSide> ret2 = snpHandler.createSNP("AACTTC", "CCGTTT", 2,
//                                                        3, halSeq,
//                                                        srcPos, tgtPos, true,
 //                                                       &lookup, NULL);
}

CuSuite* snpHandlerTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, snpMapTest);
  SUITE_ADD_TEST(suite, snpHandlerSingleSNPTest);
  SUITE_ADD_TEST(suite, snpHandlerMultibaseSNPTest);
  SUITE_ADD_TEST(suite, snpHandlerOverlapSNPTest);
  SUITE_ADD_TEST(suite, snpHandlerInversionSNPTest);
  return suite;
}
