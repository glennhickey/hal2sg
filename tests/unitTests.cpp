/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include "unitTests.h"
#include "halAlignmentTest.h"
#include "halAlignmentInstance.h"

int runAllTests(void) {
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew(); 
  CuSuiteAddSuite(suite, snpHandlerTestSuite());
  CuSuiteAddSuite(suite, sgBuildTestSuite());
  CuSuiteRun(suite);
  CuSuiteSummary(suite, output);
  CuSuiteDetails(suite, output);
  printf("%s\n", output->buffer);
  return suite->failCount > 0;
}

int main(int argc, char *argv[]) {

  storageDriverToTest = hal::STORAGE_FORMAT_HDF5;
  return runAllTests();
}
