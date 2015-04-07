/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sgbuilder.h"
#include "halBlockMapper.h"

using namespace std;
using namespace hal;

SGBuilder::SGBuilder() : _sg(0), _root(0), _lookup(0), _mapMrca(0)
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
  _lookup = NULL;
  _seqMapBack.clear();
  _mapPath.clear();
  _mapMrca = NULL;
}

const SideGraph* SGBuilder::getSideGraph() const
{
  return _sg;
}

void SGBuilder::addGenome(const Genome* genome,
                          const Sequence* sequence, 
                          hal_index_t start,
                          hal_index_t length)
{
  assert(_luMap.find(genome->getName()) == _luMap.end());

  start = std::max((hal_index_t)0, start);
  if (length <= 0)
  {
    if (sequence == NULL)
    {
      length = genome->getSequenceLength() - start;
    }
    else
    {
      length = sequence->getSequenceLength() - start;
    }
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
    const Sequence* curSequence = si->getSequence();
    for (; si != sie; si->toNext())
    {
      if (start <= curSequence->getEndPosition() &&
          (start + length-1) >= curSequence->getStartPosition())
      {
        seqNames.push_back(curSequence->getName());
      }
    }
  }
  _lookup = new SGLookup();
  _lookup->init(seqNames);
  _luMap.insert(pair<string, SGLookup*>(genome->getName(), _lookup));

  // If sequence is not NULL, start in sequence coordinates and
  // need to be converted
  if (sequence != NULL)
  {
    start += sequence->getStartPosition();
  }

  // Compute the target
  const Genome* target = getTarget(genome);
  // Update the mapping structures.  Should verify with Joel what
  // these new parameters mean. 
  if (target != NULL)
  {
    set<const Genome*> inputSet;
    inputSet.insert(genome);
    inputSet.insert(target);
    _mapMrca = getLowestCommonAncestor(inputSet);
    inputSet.clear();
    inputSet.insert(_root);
    inputSet.insert(target);
    getGenomesInSpanningTree(inputSet, _mapPath);
  }

  // Convert sequence by sequence
  for (size_t i = 0; i < seqNames.size(); ++i)
  {
    const Sequence* curSequence = genome->getSequence(seqNames[i]);
    hal_index_t endPos = start + length - 1;
    hal_index_t curStart = std::max(start, curSequence->getStartPosition());
    hal_index_t curEnd = std::min(endPos, curSequence->getEndPosition());

    // note all coordinates global
    mapSequence(curSequence, curStart, curEnd, target);  
  }
}

const Genome* SGBuilder::getTarget(const Genome* genome)
{
  // This code may be too dumb to be any more than a placeholder, but
  // have no short term designs on huge trees...

  // find the nearest genome in our hal tree that's already been added
  // to the lookup structure.
  size_t dist = numeric_limits<size_t>::max();
  const Genome* best = NULL;
  for (GenomeLUMap::const_iterator i = _luMap.begin(); i != _luMap.end(); ++i)
  {
    const Genome* candidate = _alignment->openGenome(i->first);
    size_t cdist = getGenomeDistance(candidate, genome);
    if (cdist > 0 && cdist < dist)
    {
      dist = cdist;
      best = candidate;
    }
  }
  return best;
}

SGSequence* SGBuilder::createSGSequence(const Sequence* sequence,
                                        hal_index_t startOffset,
                                        hal_index_t length)
{
  // make a new Side Graph Sequence
  SGSequence* sgSeq = new SGSequence();
  sgSeq->_length = length;
  sgSeq->_name = sequence->getName();
  // add to the Side Graph
  _sg->addSequence(sgSeq);
  // keep record of where it came from
  _seqMapBack.insert(pair<SGSequence*, pair<const Sequence*, hal_index_t> >(
                       sgSeq,
                       pair<const Sequence*, hal_index_t>(sequence,
                                                          startOffset)));
  return sgSeq;
}

void SGBuilder::mapSequence(const Sequence* sequence,
                            hal_index_t globalStart,
                            hal_index_t globalEnd,
                            const Genome* target)
{
  //easy case: add sequences without any joins because we are doing root
  //and everything is empty (same logic for empty sequence)
  const Genome* genome = sequence->getGenome();
  if ((target == NULL && genome == _root) || globalEnd == globalStart)
  {
    createSGSequence(sequence, globalStart - sequence->getStartPosition(),
                     globalEnd - globalStart + 1);
  }
  else
  {
    if (target == NULL)
    {
      // map to self
      target = genome;
      _mapMrca = genome; 
    }
    
    // block mapping logic below largely taken from halBlockLiftover.cpp
    SegmentIteratorConstPtr refSeg;
    hal_index_t lastIndex;
    if (genome->getNumTopSegments() > 0)
    {
      refSeg = genome->getTopSegmentIterator();
      lastIndex = (hal_index_t)genome->getNumTopSegments();
    }
    else
    {
      refSeg = genome->getBottomSegmentIterator();
      lastIndex = (hal_index_t)genome->getNumBottomSegments();      
    }
    
    refSeg->toSite(globalStart, false);
    hal_offset_t startOffset = globalStart - refSeg->getStartPosition();
    hal_offset_t endOffset = 0;
    if (globalEnd <= refSeg->getEndPosition())
    {
      endOffset = refSeg->getEndPosition() - globalEnd;
    }
    refSeg->slice(startOffset, endOffset);
  
    assert(refSeg->getStartPosition() ==  globalStart);
    assert(refSeg->getEndPosition() <= globalEnd);
    
    set<MappedSegmentConstPtr> mappedSegments;

    cout << "Mapping " << sequence->getFullName() << " to "
         << target->getName()
         << "; MRC = " << _mapMrca->getName()
         << "; ROOT = " << _root->getName()
         << "; PATH = ";
    for (set<const Genome*>::iterator i = _mapPath.begin();
         i != _mapPath.end(); ++i)
    {
       cout << (*i)->getName() << ", ";
    }
    cout << endl;
    
    while (refSeg->getArrayIndex() < lastIndex &&
           refSeg->getStartPosition() <= globalEnd)  
    {
      refSeg->getMappedSegments(mappedSegments, target, &_mapPath,
                                true, 0, _root, _mapMrca);
      refSeg->toRight(globalEnd);
    }

    vector<MappedSegmentConstPtr> fragments;
    BlockMapper::MSSet emptySet;
    set<hal_index_t> queryCutSet;
    set<hal_index_t> targetCutSet;

    for (set<MappedSegmentConstPtr>::iterator i = mappedSegments.begin();
         i != mappedSegments.end(); ++i)
    {
      BlockMapper::extractSegment(i, emptySet, fragments, &mappedSegments, 
                                  targetCutSet, queryCutSet);
    }
  }
}

