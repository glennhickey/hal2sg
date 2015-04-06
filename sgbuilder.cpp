/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sgbuilder.h"

using namespace std;
using namespace hal;

SGBuilder::SGBuilder() : _sg(0), _root(0)
{

}

SGBuilder::~SGBuilder()
{
  clear();
}

void SGBuilder::init(AlignmentConstPtr alignment, const Genome* root)
{
  clear();
  _alignment = alignment;
  _sg = new SideGraph();
  _root = root;
}

void SGBuilder::clear()
{
  delete _sg;
  _sg = NULL;
  _alignment = AlignmentConstPtr();
  for (GenomeLUMap::iterator i = _luMap.begin(); i != _luMap.end(); ++i)
  {
    delete i->second;
  }
  _luMap.clear();
}

void SGBuilder::addGenome(const Genome* genome,
                          const Sequence* sequence, 
                          hal_index_t start,
                          hal_index_t end)
{
  assert(_luMap.find(genome->getName()) == _luMap.end());

  if (start == NULL_INDEX)
  {
    start = 0;
  }
  if (end == NULL_INDEX)
  {
    end = genome->getSequenceLength() - 1;
  }
  
  // Add genomes sequences to the lookup structure
  vector<string> seqNames;
  if (sequence != NULL)
  {
    seqNames.push_back(sequence->getName());
  }
  else
  {
    SequenceIteratorConstPtr si = genome->getSequenceIterator();
    SequenceIteratorConstPtr sie = genome->getSequenceEndIterator();
    for (; si != sie; si->toNext())
    {
      if (start <= sequence->getEndPosition() &&
          end >= sequence->getStartPosition())
      {
        seqNames.push_back(si->getSequence()->getName());
      }
    }
  }
  SGLookup* sgLook = new SGLookup();
  sgLook->init(seqNames);
  _luMap.insert(pair<string, SGLookup*>(genome->getName(), sgLook));

  // If sequence is not NULL, start / end are sequence coordinates and
  // need to be converted
  if (sequence != NULL)
  {
    start += sequence->getStartPosition();
    end += sequence->getStartPosition();
  }

  // Compute the target
  const Genome* target = getTarget(genome);

  // Convert sequence by sequence
  for (size_t i = 0; i < seqNames.size(); ++i)
  {
    const Sequence* curSequence = genome->getSequence(seqNames[i]);
    hal_index_t curStart = std::max(start, curSequence->getStartPosition());
    hal_index_t curEnd = std::min(end, curSequence->getEndPosition());

    // note all coordinates global
    mapSequence(curSequence, curStart, curEnd, target);  
  }
}

const Genome* SGBuilder::getTarget(const Genome* genome)
{
  return NULL;
}

void SGBuilder::mapSequence(const Sequence* sequence,
                            hal_index_t globalStart,
                            hal_index_t globalEnd,
                            const Genome* target)
{

}


