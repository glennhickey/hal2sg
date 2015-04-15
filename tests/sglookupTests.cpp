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
  seqNames.push_back("seq");
  SGLookup lookup;
  lookup.init(seqNames);
  
  lookup.addInterval(SGPosition(0, 0), SGPosition(10, 10), 5, false);
  lookup.addInterval(SGPosition(0, 10), SGPosition(20, 0), 10, false);
  // try a reverse mapping
  lookup.addInterval(SGPosition(0, 20), SGPosition(50, 20), 10, true);
  
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

CuSuite* sgLookupTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, simpleTest);
  SUITE_ADD_TEST(suite, mapTest);
  return suite;
}
