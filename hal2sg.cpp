/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <deque>

#include "sgbuilder.h"
#include "sgsql.h"

using namespace std;
using namespace hal;

static bool isCamelHal(AlignmentConstPtr aligment);
static void breadthFirstGenomeSearch(const Genome* root,
                                     vector<const Genome*>& outTraversal);

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
  optionsParser->addOptionFlag("onlySequenceNames",
                               "use only sequence names "
                               "for output names.  By default, the UCSC convention of Genome.Sequence "
                               "is used",
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
  bool onlySequenceNames;
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
    onlySequenceNames = optionsParser->getFlag("onlySequenceNames");
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
    // TEMP
    if (refGenomeName != "\"\"" || targetGenomes != "\"\"" ||
        refSequenceName != "\"\"" ||
        start != 0 || length != 0 || noAncestors == true)
    {
      throw hal_exception("The following options are disabled in this release:"
                          "\n--refGenome\n--targetGenomes\n--refSeqeunce\n"
                          "--start\n--length\n--noAncestors");
    }
    
    ofstream fastaStream(fastaPath.c_str());
    if (!fastaStream)
    {
      throw hal_exception("error opening output fasta file " + fastaPath);
    }
    fastaStream.close();
    
    ofstream sqlStream(sqlPath.c_str());
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
      breadthFirstGenomeSearch(rootGenome, targetVec);
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

    bool camelMode = isCamelHal(alignment);
    if (camelMode)
    {
      cout << "CAMEL output detected.  Will infer root sequence "
           << "from children" << endl;
    }
    SGBuilder sgbuild;
    sgbuild.init(alignment, rootGenome, false, isCamelHal(alignment),
                 onlySequenceNames);
    // add the reference genome
    sgbuild.addGenome(refGenome, refSequence, start, length);

    // add the other genomes
    for (size_t i = 0; i < targetVec.size(); ++i)
    {
      sgbuild.addGenome(targetVec[i]);
    }
    
    //cout << *sgbuild.getSideGraph() << endl;

    SGSQL sqlWriter;
    sqlWriter.writeDb(&sgbuild, sqlPath, fastaPath, halPath);

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

/** CAMEL writes the root's DNA sequence as N's.  This screws up SNP detection
 * in the conversion when the root is used as the first anchor.  So we run a
 * check here to see if the root's all N's, if it is we set the camel flag
 * in sgbuilder to tell it to infer the root from its children (possible as
 * there are no substitutions in CAMEL HAL files */
bool isCamelHal(AlignmentConstPtr alignment)
{
  const Genome* rootGenome = alignment->openGenome(alignment->getRootName());
  if (rootGenome->getSequenceLength() == 0)
  {
    return false;
  }
  DNAIteratorConstPtr dnaIt = rootGenome->getDNAIterator();
  DNAIteratorConstPtr end = rootGenome->getDNAEndIterator();
  for (; !dnaIt->equals(end); dnaIt->toRight())
  {
    if (toupper(dnaIt->getChar()) != 'N')
    {
      return false;
    }
  }
  return true;
}

void breadthFirstGenomeSearch(const Genome* root,
                              vector<const Genome*>& outTraversal)
{
  deque<const Genome*> bfsQueue;
  bfsQueue.push_back(root);
  while (!bfsQueue.empty())
  {
    const Genome* genome = bfsQueue.front();
    bfsQueue.pop_front();
    outTraversal.push_back(genome);
    for (hal_size_t i = 0; i < genome->getNumChildren(); ++i)
    {
      bfsQueue.push_back(genome->getChild(i));
    }
  }
}
