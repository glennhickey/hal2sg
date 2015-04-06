/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>

#include "hal.h"

using namespace std;
using namespace hal;



static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("sgFile", "SQL info here todo..");
  optionsParser->addOption("refGenome", 
                           "name of reference genome (root if empty)", 
                           "\"\"");
  optionsParser->addOption("refSequence",
                           "name of reference sequence within reference genome"
                           " (all sequences if empty)",
                           "\"\"");
  optionsParser->addOption("start",
                           "coordinate within reference genome (or sequence"
                           " if specified) to start at",
                           0);
  optionsParser->addOption("length",
                           "length of the reference genome (or sequence"
                           " if specified) to convert.  If set to 0,"
                           " the entire thing is converted",
                           0);
  optionsParser->addOption("rootGenome", 
                           "name of root genome (none if empty)", 
                           "\"\"");
  optionsParser->addOption("targetGenomes",
                           "comma-separated (no spaces) list of target genomes "
                           "(others are excluded) (vist all if empty)",
                           "\"\"");
  optionsParser->addOptionFlag("noAncestors", 
                               "don't write ancestral sequences. IMPORTANT: "
                               "Must be used in conjunction with --refGenome"
                               " to set a non-ancestral genome as the reference"
                               " because the default reference is the root.", 
                               false);

  optionsParser->setDescription("Convert hal database to GA4GH Side Graph");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();
  string halPath;
  string refGenomeName;
  string rootGenomeName;
  string targetGenomes;
  string refSequenceName;
  hal_index_t start;
  hal_size_t length;
  bool noAncestors;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    refGenomeName = optionsParser->getOption<string>("refGenome");
    rootGenomeName = optionsParser->getOption<string>("rootGenome");
    targetGenomes = optionsParser->getOption<string>("targetGenomes");
    refSequenceName = optionsParser->getOption<string>("refSequence");    
    start = optionsParser->getOption<hal_index_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
    noAncestors = optionsParser->getFlag("noAncestors");

    if (rootGenomeName != "\"\"" && targetGenomes != "\"\"")
    {
      throw hal_exception("--rootGenome and --targetGenomes options are "
                          "mutually exclusive");
    }
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }
  try
  {
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath, 
                                                           optionsParser);
    if (alignment->getNumGenomes() == 0)
    {
      throw hal_exception("hal alignmenet is empty");
    }
    
    set<const Genome*> targetSet;
    const Genome* rootGenome = NULL;
    if (rootGenomeName != "\"\"")
    {
      rootGenome = alignment->openGenome(rootGenomeName);
      if (rootGenome == NULL)
      {
        throw hal_exception(string("Root genome, ") + rootGenomeName + 
                            ", not found in alignment");
      }
      if (rootGenomeName != alignment->getRootName())
      {
        getGenomesInSubTree(rootGenome, targetSet);
      }
    }

    if (targetGenomes != "\"\"")
    {
      vector<string> targetNames = chopString(targetGenomes, ",");
      for (size_t i = 0; i < targetNames.size(); ++i)
      {
        const Genome* tgtGenome = alignment->openGenome(targetNames[i]);
        if (tgtGenome == NULL)
        {
          throw hal_exception(string("Target genome, ") + targetNames[i] + 
                              ", not found in alignment");
        }
        targetSet.insert(tgtGenome);
      }
    }

    const Genome* refGenome = NULL;
    if (refGenomeName != "\"\"")
    {
      refGenome = alignment->openGenome(refGenomeName);
      if (refGenome == NULL)
      {
        throw hal_exception(string("Reference genome, ") + refGenomeName + 
                            ", not found in alignment");
      }
    }
    else
    {
      refGenome = alignment->openGenome(alignment->getRootName());
    }
    const SegmentedSequence* ref = refGenome;

    if (noAncestors == true && refGenome->getNumChildren() != 0)
    {
      throw hal_exception(string("Since the reference genome to be used for the"
                                 " MAF is ancestral (") + refGenome->getName() +
                          "), the --noAncestors option is invalid.  The "
                          "--refGenome option can be used to specify a "
                          "different reference.");
    }
    
    const Sequence* refSequence = NULL;
    if (refSequenceName != "\"\"")
    {
      refSequence = refGenome->getSequence(refSequenceName);
      ref = refSequence;
      if (refSequence == NULL)
      {
        throw hal_exception(string("Reference sequence, ") + refSequenceName + 
                            ", not found in reference genome, " + 
                            refGenome->getName());
      }
    }
  }
  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
  
  return 0;
}
