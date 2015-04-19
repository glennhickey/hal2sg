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

  SNPHandler snpHandler(false);

  snpHandler.addSNP(p1, 'A', p4);
  snpHandler.addSNP(p3, 'A', p2);
  snpHandler.addSNP(p1, 'c', p3);
  snpHandler.addSNP(p1, 'T', p2);

  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'A') == p4);
  CuAssertTrue(testCase, snpHandler.findSNP(p3, 'A') == p2);
  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'c') == p3);
  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'C') == p3);
  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'T') == p2);
  CuAssertTrue(testCase, snpHandler.findSNP(p1, 'G') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandler.findSNP(p3, 'G') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandler.findSNP(p2, 'A') == SideGraph::NullPos);

  SNPHandler snpHandlerCS(true);

  snpHandlerCS.addSNP(p1, 'A', p4);
  snpHandlerCS.addSNP(p3, 'A', p2);
  snpHandlerCS.addSNP(p1, 'c', p3);
  snpHandlerCS.addSNP(p1, 'T', p2);
  snpHandlerCS.addSNP(p1, 't', p3);

  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'A') == p4);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p3, 'A') == p2);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'c') == p3);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'C') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'T') == p2);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 't') == p3);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p1, 'G') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p3, 'G') == SideGraph::NullPos);
  CuAssertTrue(testCase, snpHandlerCS.findSNP(p2, 'A') == SideGraph::NullPos);
}

CuSuite* snpHandlerTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, snpMapTest);
  return suite;
}
