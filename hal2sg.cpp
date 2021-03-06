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
#include "halsgsql.h"

using namespace std;
using namespace hal;

static bool isCamelHal(AlignmentConstPtr aligment);
static void breadthFirstGenomeSearch(const Genome* reference,
                                     const vector<const Genome*>& targets,
                                     vector<const Genome*>& outTraversal);


static void initParser(CLParser* optionsParser)
{
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("fastaFile", "Output FASTA sequences");
  optionsParser->addArgument("sqlFile", "SQL inserts written here");
  optionsParser->addOption("refGenome", 
                           "name of reference genome (HAL root if empty)", 
                           "\"\"");
  optionsParser->addOption("rootGenome", 
                           "process only genomes in clade with specified root"
                           " (HAL root if empty)", 
                           "\"\"");
  optionsParser->addOption("targetGenomes",
                           "comma-separated (no spaces) list of target genomes "
                           "(others are excluded) (vist all if empty)",
                           "\"\"");
  optionsParser->addOptionFlag("noAncestors", 
                               "don't write ancestral paths. Note that "
                               "ancestral *sequences* may still get written"
                               " as they can be necessary for expressing"
                               " some alignments. IMPORTANT: "
                               "Must be used in conjunction with --refGenome"
                               " to set a non-ancestral genome as the reference"
                               " because the default reference is the root.", 
                               false);
  optionsParser->addOptionFlag("onlySequenceNames",
                               "use only sequence names for output names.  By "
                               "default, the UCSC convention of "
                               "Genome.Sequence is used",
                               false);

  optionsParser->setDescription("Convert HAL alignment to GA4GH Side "
                                "Graph SQL format");
}

int main(int argc, char** argv)
{
  CLParser optionsParser;
  initParser(&optionsParser);
  string halPath;
  string fastaPath;
  string sqlPath;
  string refGenomeName;
  string rootGenomeName;
  string targetGenomes;
  bool noAncestors;
  bool onlySequenceNames;
  try
  {
    optionsParser.parseOptions(argc, argv);
    halPath = optionsParser.getArgument<string>("halFile");
    fastaPath = optionsParser.getArgument<string>("fastaFile");
    sqlPath = optionsParser.getArgument<string>("sqlFile");
    refGenomeName = optionsParser.getOption<string>("refGenome");
    rootGenomeName = optionsParser.getOption<string>("rootGenome");
    targetGenomes = optionsParser.getOption<string>("targetGenomes");
    noAncestors = optionsParser.getFlag("noAncestors");
    onlySequenceNames = optionsParser.getFlag("onlySequenceNames");
    if (rootGenomeName != "\"\"" && targetGenomes != "\"\"")
    {
      throw hal_exception("--rootGenome and --targetGenomes options are "
                          "mutually exclusive");
    }
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser.printUsage(cerr);
    exit(1);
  }
  try
  {
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
    
    AlignmentConstPtr alignment(openHalAlignment(halPath, 
                                                 &optionsParser,
                                                 hal::READ_ACCESS));
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

    // if target set not specified we default to all leaves under the
    // given root. 
    vector<const Genome*> targetVec;
    if (targetGenomes == "\"\"")
    {
      set<const Genome*> allGenomes;
      getGenomesInSubTree(rootGenome, allGenomes);
      for (set<const Genome*>::iterator i = allGenomes.begin();
           i != allGenomes.end(); ++i)
      {
        if ((*i)->getNumChildren() == 0)
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

    // get a breadth-first ordering of all target genomes
    // starting at the reference and including any ancestors
    // that need to get walked across
    vector<const Genome*> breadthFirstOrdering;
    breadthFirstGenomeSearch(refGenome, targetVec, breadthFirstOrdering);

    bool camelMode = isCamelHal(alignment);
    if (camelMode)
    {
      cout << "CAMEL output detected.  Will infer root sequence "
           << "from children" << endl;
    }
    SGBuilder sgbuild;
    sgbuild.init(alignment, rootGenome, false, isCamelHal(alignment),
                 onlySequenceNames);
    
    // add the genomes in the breadth first order
    for (size_t i = 0; i < breadthFirstOrdering.size(); ++i)
    {
      sgbuild.addGenome(breadthFirstOrdering[i]);
    }

    // compute all the joins in second pass (and do sanity check
    // on every path in graph)
    sgbuild.computeJoins(!noAncestors);
    
    //cout << *sgbuild.getSideGraph() << endl;

    HALSGSQL sqlWriter;
    sqlWriter.exportGraph(&sgbuild, sqlPath, fastaPath, halPath,
                          !noAncestors);

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
  DnaIteratorPtr dnaIt = rootGenome->getDnaIterator();
  size_t length = rootGenome->getSequenceLength();
  for (hal_size_t i = 0; i < length; ++i, dnaIt->toRight())
  {
    if (toupper(dnaIt->getBase()) != 'N')
    {
      return false;
    }
  }
  return true;
}

void breadthFirstGenomeSearch(const Genome* reference,
                              const vector<const Genome*>& targets,
                              vector<const Genome*>& outTraversal)
{
  // find all genomes we need to visit
  set<const Genome*> inputSet;
  inputSet.insert(reference);
  for (size_t i = 0; i < targets.size(); ++i)
  {
    inputSet.insert(targets[i]);
  }
  set<const Genome*> visitSet;
  getGenomesInSpanningTree(inputSet, visitSet);

  // find our breadth first order through the visit set, starting at
  // reference.
  set<const Genome*> flagged;
  deque<const Genome*> bfsQueue;
  bfsQueue.push_back(reference);
  while (!bfsQueue.empty())
  {
    const Genome* genome = bfsQueue.front();
    bfsQueue.pop_front();
    outTraversal.push_back(genome);
    vector<const Genome*> neighbours;
    for (hal_size_t i = 0; i < genome->getNumChildren(); ++i)
    {
      neighbours.push_back(genome->getChild(i));
    }
    if (genome->getParent() != NULL)
    {
      neighbours.push_back(genome->getParent());
    }
    for (hal_size_t i = 0; i < neighbours.size(); ++i)
    {
      const Genome* neighbour = neighbours[i];
      if (visitSet.find(neighbour) != visitSet.end() &&
          flagged.find(neighbour) == flagged.end())
      {
        bfsQueue.push_back(neighbour);
      }
    }
    flagged.insert(genome);
  }
}
