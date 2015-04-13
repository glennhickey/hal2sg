/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>

#include "sgbuilder.h"
#include "sgsql.h"

using namespace std;
using namespace hal;



static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("fastaFile", "Output FASTA sequences");
  optionsParser->addArgument("sqlFile", "SQL inserts written here");
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
  string fastaPath;
  string sqlPath;
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
    fastaPath = optionsParser->getArgument<string>("fastaFile");
    sqlPath = optionsParser->getArgument<string>("sqlFile");
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
    ofstream fastaStream(fastaPath);
    if (!fastaStream)
    {
      throw hal_exception("error opening output fasta file " + fastaPath);
    }
    fastaStream.close();
    
    ofstream sqlStream(sqlPath);
    if (!sqlStream)
    {
      throw hal_exception("error opening output sql file " + sqlPath);
    }
    sqlStream.close();
    
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath, 
                                                           optionsParser);
    if (alignment->getNumGenomes() == 0)
    {
      throw hal_exception("hal alignmenet is empty");
    }

    // root is specified either by the parameter or as the alignment root
    // by default
    const Genome* rootGenome = NULL;
    if (rootGenomeName != "\"\"")
    {
      rootGenome = alignment->openGenome(rootGenomeName);
    }
    else
    {
      rootGenome = alignment->openGenome(alignment->getRootName());
    }
    if (rootGenome == NULL)
    {
      throw hal_exception(string("Root genome, ") + rootGenomeName + 
                          ", not found in alignment");
    }

    // target genomes pulled from tree traversal (using optional root
    // parameter)
    vector<const Genome*> targetVec;
    if (targetGenomes == "\"\"")
    {
      set<const Genome*> targetSet;
      getGenomesInSubTree(rootGenome, targetSet);
      for (set<const Genome*>::iterator i = targetSet.begin();
           i != targetSet.end(); ++i)
      {
        if (noAncestors == false || (*i)->getNumChildren() == 0)
        {
          targetVec.push_back(*i);
        }
      }
    }
    // target genomes pulled from list.  
    else
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
        targetVec.push_back(tgtGenome);
      }
      
    }

    // open the reference genome (root genome if unspecified)
    const Genome* refGenome = NULL;
    if (refGenomeName != "\"\"")
    {
      refGenome = alignment->openGenome(refGenomeName);
      if (refGenome == NULL)
      {
        throw hal_exception(string("Reference genome, ") + refGenomeName + 
                            ", not found in alignment");
      }
      set<const Genome*> genomeSet;
      genomeSet.insert(refGenome);
      genomeSet.insert(rootGenome);
      if (getLowestCommonAncestor(genomeSet) != rootGenome)
      {
        throw hal_exception(string("reference genome must be under root"));
      }
    }
    else
    {
      refGenome = rootGenome;
    }

    // optionally specify a sequence in the ref genome
    const Sequence* refSequence = NULL;
    if (refSequenceName != "\"\"")
    {
      refSequence = refGenome->getSequence(refSequenceName);
      if (refSequence == NULL)
      {
        throw hal_exception(string("Reference sequence, ") + refSequenceName + 
                            ", not found in reference genome, " + 
                            refGenome->getName());
      }
    }

    // make sure refGenome not in target genomes
    for (vector<const Genome*>::iterator i = targetVec.begin();
         i != targetVec.end(); ++i)
    {
      if (*i == refGenome)
      {
        targetVec.erase(i);
        break;
      }
    }

    SGBuilder sgbuild;
    sgbuild.init(alignment, rootGenome, false, true);
    // add the reference genome
    sgbuild.addGenome(refGenome, refSequence, start, length);

    // add the other genomes
    for (size_t i = 0; i < targetVec.size(); ++i)
    {
      sgbuild.addGenome(targetVec[i]);
    }
    
    cout << *sgbuild.getSideGraph() << endl;

    SGSQL sqlWriter;
    sqlWriter.writeDb(&sgbuild, sqlPath, fastaPath);

  }
/*  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
*/
  catch(int e) {}
  return 0;
}
