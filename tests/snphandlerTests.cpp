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
                                                    false, false, &lookup,
                                                    NULL);

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
                                                    false, false,
                                                    &lookup, NULL);

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
  string dna1 = "T123ATA7890TT";
  string dna2 = "TC....AGGGGGA";
  string dnaTgt = "GGGGGGGTTTTAG";

  cout << "\n ROW 1" << endl;
  SGPosition srcPos(1, 0);
  SGPosition tgtPos(0, 0);
  snpHandler.createSNP(dna0, dnaTgt, 0, 2, halSeq, srcPos, tgtPos, false,
                       false, &lookup, NULL);
  srcPos.setPos(5);
  tgtPos.setPos(5);
  snpHandler.createSNP(dna0, dnaTgt, 5, 1, halSeq, srcPos, tgtPos, false,
                       false, &lookup, NULL);
  srcPos.setPos(9);
  tgtPos.setPos(9);
  snpHandler.createSNP(dna0, dnaTgt, 9, 3, halSeq, srcPos, tgtPos, false,
                       false, &lookup, NULL);

  cout << "\n ROW 2" << endl;
  srcPos = SGPosition(2, 0);
  tgtPos = SGPosition(0, 0);
  snpHandler.createSNP(dna1, dnaTgt, 0, 1, halSeq, srcPos, tgtPos, false,
                       false, &lookup, NULL);
  srcPos.setPos(4);
  tgtPos.setPos(4);
  snpHandler.createSNP(dna1, dnaTgt, 4, 3, halSeq, srcPos, tgtPos, false,
                       false, &lookup, NULL);
  srcPos.setPos(11);
  tgtPos.setPos(11);
  snpHandler.createSNP(dna1, dnaTgt, 11, 2, halSeq, srcPos, tgtPos, false,
                       false, &lookup, NULL);

  cout << "\n ROW 3" << endl;
  srcPos = SGPosition(3, 0);
  tgtPos = SGPosition(0, 0);
  snpHandler.createSNP(dna2, dnaTgt, 0, 2, halSeq, srcPos, tgtPos, false,
                       false, &lookup, NULL);
  srcPos.setPos(6);
  tgtPos.setPos(6);
  snpHandler.createSNP(dna2, dnaTgt, 6, 7, halSeq, srcPos, tgtPos, false,
                       false, &lookup, NULL);

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
  snpHandler.createSNP("A", "G", 0, 1, halSeq,
                       srcPos, tgtPos, true,
                       false,
                       &lookup, NULL);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,50), 'G').getSeqID() == 0);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,50), 'T').getSeqID() == 1);

  srcPos.setPos(150);
  tgtPos.setPos(150);
  snpHandler.createSNP("A", "G", 0, 1, halSeq,
                       srcPos, tgtPos, false,
                       true,
                       &lookup, NULL);
  // note: unlike above, our baseline is reversed because tgt value is
  // reverse mapped to the side graph (sgReverseMap true).  this doesn't
  // bear on the src value, as it's still reversed (tranMap true)
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,150), 'C').getSeqID() == 0);
  CuAssertTrue(tc, snpHandler.findSNP(SGPosition(0,150), 'T').getSeqID() == 2);


  cout << "\n\n mult test \n";
  // multi SNP reverse map
  srcPos.setPos(60);
  tgtPos.setPos(60);
  snpHandler.createSNP("AACTTC", "CCGTTT", 2,
                       3, halSeq,
                       srcPos, tgtPos, true,
                       false,
                       &lookup, NULL);
  // check that we added a 3-length snp of rev(SRC) along the
  // side graph
  SGPosition s = snpHandler.findSNP(SGPosition(1,60), 'C');
  CuAssertTrue(tc, s == SideGraph::NullPos);
  s = snpHandler.findSNP(SGPosition(0, 60), 'G');
  CuAssertTrue(tc, s.getPos() == 2);
  CuAssertTrue(tc, s.getSeqID() == 3);
  s = snpHandler.findSNP(SGPosition(0, 59), 'A');
  CuAssertTrue(tc, s.getPos() == 1);
  CuAssertTrue(tc, s.getSeqID() == 3);
  s = snpHandler.findSNP(SGPosition(0, 58), 'A');
  CuAssertTrue(tc, s.getPos() == 0);
  CuAssertTrue(tc, s.getSeqID() == 3);

  cout << endl;
  // check that we added the baseline snps of the side graph
  s = snpHandler.findSNP(SGPosition(0, 58), 'C');
  CuAssertTrue(tc, s.getPos() == 58);
  CuAssertTrue(tc, s.getSeqID() == 0);
  s = snpHandler.findSNP(SGPosition(0, 59), 'G');
  CuAssertTrue(tc, s.getPos() == 59);
  CuAssertTrue(tc, s.getSeqID() == 0);
  s = snpHandler.findSNP(SGPosition(0, 60), 'T');
  CuAssertTrue(tc, s.getPos() == 60);
  CuAssertTrue(tc, s.getSeqID() == 0);

  // case inspired by problem mapping ref:794727-8 to Anc0 in sma_star.hal
  srcPos.setPos(1000);
  tgtPos.setPos(1000);
  snpHandler.createSNP("CCCAT", "AATGG",
                       2, 2, halSeq, srcPos, tgtPos, true, false,
                       &lookup, NULL);
  s = snpHandler.findSNP(SGPosition(0,1000), 'G');
  CuAssertTrue(tc, s.getPos() == 1);
  CuAssertTrue(tc, s.getSeqID() == 4);
  s = snpHandler.findSNP(SGPosition(0, 999), 'T');
  CuAssertTrue(tc, s.getPos() == 0);
  CuAssertTrue(tc, s.getSeqID() == 4);


  // another case
  srcPos.setPos(1500);
  tgtPos.setPos(1500);
  snpHandler.createSNP("T", "C", 0, 1, halSeq,
                       SGPosition(srcPos.getSeqID(), 2500),
                       SGPosition(tgtPos.getSeqID(), 1500),
                       false, false, &lookup, NULL);

  snpHandler.createSNP("AC", "AC", 0, 2, halSeq, srcPos, tgtPos, true, false,
                       &lookup, NULL);
                     
  s = snpHandler.findSNP(SGPosition(0,1500), 'T');
  CuAssertTrue(tc, s.getPos() == 0);
  CuAssertTrue(tc, s.getSeqID() == 5);
  s = snpHandler.findSNP(SGPosition(0, 1499), 'G');
  CuAssertTrue(tc, s.getPos() == 0);
  CuAssertTrue(tc, s.getSeqID() == 6);

  // same as above, but with substrings case
  srcPos.setPos(2000);
  tgtPos.setPos(2000);
  string sdna(53, '.');
  string tdna(53, ',');
  sdna[40] = 'A';
  sdna[41] = 'C';
  tdna[53 - 1 - 41] = 'A';
  tdna[53 - 1 - 40] = 'C';
  snpHandler.createSNP("T", "A", 0, 1, halSeq,
                       SGPosition(srcPos.getSeqID(), 5),
                       SGPosition(tgtPos.getSeqID(), 2000),
                       false, false, &lookup, NULL);
                     
  snpHandler.createSNP(sdna, tdna, 40, 2, halSeq, srcPos, tgtPos, true, false,
     &lookup, NULL);
  s = snpHandler.findSNP(SGPosition(0,2000), 'T');
  CuAssertTrue(tc, s.getPos() == 0);
  CuAssertTrue(tc, s.getSeqID() == 7);
  s = snpHandler.findSNP(SGPosition(0, 1999), 'G');
  CuAssertTrue(tc, s.getPos() == 0);
  CuAssertTrue(tc, s.getSeqID() == 8);

  // similar deal but with multiple snps already added
  srcPos.setPos(5000);
  tgtPos.setPos(5000);

  string dna1 = "GCAAATGGCGGTC";
  string dna2 = "GAATAGTTAGTGC";

  cout << "\n\n";
  snpHandler.createSNP(dna1, dna2, 4, 2, halSeq,
                       SGPosition(srcPos.getSeqID(), 3004),
                       SGPosition(tgtPos.getSeqID(), 4996),
                       true, false, &lookup, NULL);

  snpHandler.createSNP(dna1, dna2, 8, 1, halSeq,
                       SGPosition(srcPos.getSeqID(), 3008),
                       SGPosition(tgtPos.getSeqID(), 4992),
                       true, false, &lookup, NULL);

  snpHandler.createSNP(dna1, dna2, 3, 8, halSeq,
                       SGPosition(srcPos.getSeqID(), 5003),
                       SGPosition(tgtPos.getSeqID(), 4997),
                       true, false, &lookup, NULL);


  
  s = snpHandler.findSNP(SGPosition(0,1500), 'T');
  CuAssertTrue(tc, s.getPos() == 0);
  CuAssertTrue(tc, s.getSeqID() == 5);
  s = snpHandler.findSNP(SGPosition(0, 1499), 'G');
  CuAssertTrue(tc, s.getPos() == 0);
  CuAssertTrue(tc, s.getSeqID() == 6);
  

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
