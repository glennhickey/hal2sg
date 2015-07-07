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
#include "sglookup.h"

using namespace std;

/** randome test below too hard to start.. */
void simpleTest(CuTest *testCase)
{
  vector<string> seqNames;
  seqNames.push_back("seq1");
  seqNames.push_back("seq2");
  SGLookup lookup;
  lookup.init(seqNames);
  
  lookup.addInterval(SGPosition(0, 0), SGPosition(10, 10), 5, false);
  lookup.addInterval(SGPosition(0, 10), SGPosition(20, 0), 10, false);
  // try a reverse mapping
  lookup.addInterval(SGPosition(0, 20), SGPosition(50, 20), 10, true);

  // try a reverse mapping followed by forward
  lookup.addInterval(SGPosition(1, 0), SGPosition(60, 0), 10, true);
  lookup.addInterval(SGPosition(1, 10), SGPosition(60, 10), 30, false);
  lookup.addInterval(SGPosition(1, 40), SGPosition(60, 40), 5, false);
  
  for (size_t i = 0; i < 5; ++i)
  {
    SGPosition mapPos = lookup.mapPosition(SGPosition(0, i)).getBase();
    CuAssertTrue(testCase, mapPos == SGPosition(10, 10 + i));
  }
  for (size_t i = 0; i < 10; ++i)
  {
    SGSide mapSide = lookup.mapPosition(SGPosition(0, 10 + i));
    SGPosition mapPos = mapSide.getBase();
    CuAssertTrue(testCase, mapSide.getForward() == true);
    CuAssertTrue(testCase, mapPos == SGPosition(20, i));

  }
  for (size_t i = 5; i < 10; ++i)
  {
    SGSide mapSide = lookup.mapPosition(SGPosition(0, i));
    SGPosition mapPos = mapSide.getBase();
    CuAssertTrue(testCase, mapSide.getForward() == true);
    CuAssertTrue(testCase, mapPos == SGPosition(-1, -1));
  }
  for (size_t i = 0; i < 10; ++i)
  {
    SGSide mapSide = lookup.mapPosition(SGPosition(0, 20 + i));
    SGPosition mapPos = mapSide.getBase();
    CuAssertTrue(testCase, mapSide.getForward() == false);
    CuAssertTrue(testCase, mapPos == SGPosition(50, 29-i));
  }

  for (size_t i = 0; i < 10; ++i)
  {
    SGSide mapSide = lookup.mapPosition(SGPosition(1, i));
    SGPosition mapPos = mapSide.getBase();
    CuAssertTrue(testCase, mapSide.getForward() == false);
    CuAssertTrue(testCase, mapPos == SGPosition(60, 9-i));
  }

  for (size_t i = 0; i < 35; ++i)
  {
    SGSide mapSide = lookup.mapPosition(SGPosition(1, 10 + i));
    SGPosition mapPos = mapSide.getBase();
    CuAssertTrue(testCase, mapSide.getForward() == true);
    CuAssertTrue(testCase, mapPos == SGPosition(60, 10 + i));
  }
}

/** test underlying mapping structure in lookup with a bunch of 
 * random intervals.  note: only forward strand map done for now */
void mapTest(CuTest *testCase)
{
  //srand(time(NULL));

  size_t numSequences = 50;
  size_t maxSeqLen = 100000;
  size_t maxMapSize = 1000;
  
  vector<string> seqNames;
  for (size_t i = 0; i < numSequences; ++i)
  {
    stringstream ss;
    ss << "seq" << i;
    seqNames.push_back(ss.str());
  }

  SGLookup lookup;
  lookup.init(seqNames);

  vector<vector<SGPosition> > truth;
  for (size_t i = 0; i < numSequences; ++i)
  {
    // make a source sequence
    size_t seqLen = 1 + rand() % maxSeqLen;
    truth.push_back(vector<SGPosition>(seqLen + maxMapSize + 1));
    size_t j = 0;
    while (j < seqLen)
    {
      // choose source range
      size_t k = j + 1 + rand() % maxMapSize;

      // leave holes sometimes
      if (j % 2 != 0 && i % 2 != 0 && i != 0)
      {
        // invent random target (since we do no checks on this)
        sg_seqid_t tgtID = rand() % 1000;
        sg_int_t tgtOffset = rand() % 1000000;
      
        // add block to map
        SGPosition srcPos(i, j);
        SGPosition tgtPos(tgtID, tgtOffset);
        sg_int_t length = k - j + 1;
        lookup.addInterval(srcPos, tgtPos, length, false);

        // add the block to truth
        for (sg_int_t t = 0; t < length; ++t)
        {
          truth[i][j + t] = SGPosition(tgtID, tgtOffset + t);
        }

      }

      j = k + 1;
    }
  }

  // check our lookup vs the truth
  for (size_t i = 0; i < numSequences; ++i)
  {
    for (size_t j = 0; j < truth[i].size(); ++j)
    {
      SGPosition srcPos(i, j);
      SGPosition mapPos = lookup.mapPosition(srcPos).getBase();
      CuAssertTrue(testCase, mapPos == truth[i][j]);
    }
  }
}


/** 
 * A simple path test.  Note that path code will be better tested by
 * the consistency check in sgbuilder (verifyPath())
 */
void pathTest(CuTest *testCase)
{
  vector<string> seqNames;
  seqNames.push_back("seq");
  SGLookup lookup;
  lookup.init(seqNames);
  
  vector<SGSegment> srcSegments;
  vector<SGSegment> tgtSegments;

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 0), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(20, 5), true), 10));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 10), true), 1));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(20, 4), true), 1));
  
  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 11), true), 5));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(21, 0), true), 5));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 16), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(21, 10), false), 10));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 26), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(22, 31), false), 10));
  
  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 36), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(22, 20), false), 10));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 46), true), 7));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(22, 10), true), 7));
  
  for (size_t i = 0; i < srcSegments.size(); ++i)
  {
    SGPosition halPos = srcSegments[i].getSide().getBase();
    sg_int_t length = srcSegments[i].getLength();
    SGPosition leftTgtPos = tgtSegments[i].getMinPos();
    bool reversed = tgtSegments[i].getSide().getForward() == false;
    lookup.addInterval(halPos, leftTgtPos, length, reversed);
  }

  vector<SGSegment> path;
  SGSegment last = srcSegments.back();
  last.flip();
  lookup.getPath(srcSegments[0].getSide().getBase(), last.getSide().getBase(),
                 path);

  CuAssertTrue(testCase, path == tgtSegments);

  
  // try the same, but reverse all mappings
  lookup.init(seqNames);
  srcSegments.clear();
  tgtSegments.clear();
  
  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 0), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(20, 15), false), 10));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 10), true), 1));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(20, 14), false), 1));
  
  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 11), true), 5));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(21, 10), false), 5));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 16), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(21, 110), true), 10));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 26), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(22, 130), true), 10));
  
  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 36), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(22, 120), true), 10));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 46), true), 7));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(22, 110), false), 7));

  for (size_t i = 0; i < srcSegments.size(); ++i)
  {
    SGPosition halPos = srcSegments[i].getSide().getBase();
    sg_int_t length = srcSegments[i].getLength();
    SGPosition leftTgtPos = tgtSegments[i].getMinPos();
    bool reversed = tgtSegments[i].getSide().getForward() == false;
    lookup.addInterval(halPos, leftTgtPos, length, reversed);
  }


  path.clear();
  last = srcSegments.back();
  last.flip();
  lookup.getPath(srcSegments[0].getSide().getBase(), last.getSide().getBase(),
                 path);

  CuAssertTrue(testCase, path == tgtSegments);

  // try a couple subpaths.

  path.clear();
  lookup.getPath(SGPosition(0, 5), SGPosition(0, 15), path);
  CuAssertTrue(testCase, path.size() == 3);
  CuAssertTrue(testCase, path[0] ==
               SGSegment(SGSide(SGPosition(20, 10), false), 5));
  CuAssertTrue(testCase, path[1] ==
               SGSegment(SGSide(SGPosition(20, 14), false), 1));
  CuAssertTrue(testCase, path[2] ==
               SGSegment(SGSide(SGPosition(21, 10), false), 5));


  path.clear();
  lookup.getPath(SGPosition(0, 37), SGPosition(0, 27), path);
  CuAssertTrue(testCase, path.size() == 2);
  CuAssertTrue(testCase, path[0] ==
               SGSegment(SGSide(SGPosition(22, 121), false), 2));
  CuAssertTrue(testCase, path[1] == 
               SGSegment(SGSide(SGPosition(22, 139), false), 9));

}

/** 
 * Test the outDist parameter added to mapInterval function.. 
 */
void outDistTest(CuTest *testCase)
{
  vector<string> seqNames;
  seqNames.push_back("seq");
  SGLookup lookup;
  lookup.init(seqNames);
  
  vector<SGSegment> srcSegments;
  vector<SGSegment> tgtSegments;

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 0), true), 10));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(100, 10), true), 10));

  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 10), true), 7));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(100, 50), false), 7));
  
  srcSegments.push_back(SGSegment(SGSide(SGPosition(0, 17), true), 15));
  tgtSegments.push_back(SGSegment(SGSide(SGPosition(100, 80), true), 15));

  for (size_t i = 0; i < srcSegments.size(); ++i)
  {
    SGPosition halPos = srcSegments[i].getSide().getBase();
    sg_int_t length = srcSegments[i].getLength();
    SGPosition leftTgtPos = tgtSegments[i].getMinPos();
    bool reversed = tgtSegments[i].getSide().getForward() == false;
    lookup.addInterval(halPos, leftTgtPos, length, reversed);
  }

  for (size_t i = 0; i < srcSegments.size(); ++i)
  {
    sg_int_t length = srcSegments[i].getLength();
    for (size_t j = 0; j < length; ++j)
    {
      sg_int_t dist = -1;
      SGPosition halPos = srcSegments[i].getSide().getBase();
      halPos.setPos(halPos.getPos() + j);
      SGSide result = lookup.mapPosition(halPos, &dist, false);
      CuAssertTrue(testCase, dist == length - j - 1);

      result = lookup.mapPosition(halPos, &dist, true);
      CuAssertTrue(testCase, dist == j);
    }    
  }
}

/** 
 * Add this test to help debug crash on what seems to be a fairly 
 * straightwoard case...
 */
void path2Test(CuTest *testCase)
{
  vector<string> seqNames;
  seqNames.push_back("seq0");
  seqNames.push_back("seq1");
  seqNames.push_back("seq2");
  seqNames.push_back("seq3");
  SGLookup lookup;
  lookup.init(seqNames);

  lookup.addInterval(SGPosition(0, 0), SGPosition(0,63130), 16458, 1);
  vector<SGSegment> path;
  lookup.getPath(SGPosition(0, 0), SGPosition(0, 16457), path);

  // problem seems to be bug in case where path maps to single segment
  // in reverse strand.  now fixed.
  
  CuAssertTrue(testCase, path.size() == 1);
  CuAssertTrue(testCase, path[0].getLength() == 16458);
  CuAssertTrue(testCase, path[0].getSide().getBase() ==
               SGPosition(0, 63130 + 16458 - 1));
  CuAssertTrue(testCase, path[0].getSide().getForward() == false);

  lookup.addInterval(SGPosition(1, 0), SGPosition(1,63130), 16458, 0);
  path.clear();
  lookup.getPath(SGPosition(1, 0), SGPosition(1, 16457), path);

  CuAssertTrue(testCase, path.size() == 1);
  CuAssertTrue(testCase, path[0].getLength() == 16458);
  CuAssertTrue(testCase, path[0].getSide().getBase() ==
               SGPosition(1, 63130));
  CuAssertTrue(testCase, path[0].getSide().getForward() == true);  
}

CuSuite* sgLookupTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, simpleTest);
  SUITE_ADD_TEST(suite, mapTest);
  SUITE_ADD_TEST(suite, pathTest);
  SUITE_ADD_TEST(suite, outDistTest);
  SUITE_ADD_TEST(suite, path2Test);
  return suite;
}
