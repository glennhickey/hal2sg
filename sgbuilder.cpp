/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <sstream>

#include "sgbuilder.h"
#include "halBlockMapper.h"

using namespace std;
using namespace hal;

SGBuilder::SGBuilder() : _sg(0), _root(0), _lookup(0), _mapMrca(0),
                         _referenceDupes(true)
{

}

SGBuilder::~SGBuilder()
{
  clear();
}

void SGBuilder::init(AlignmentConstPtr alignment, const Genome* root,
                     bool referenceDupes)
{
  clear();
  _alignment = alignment;
  _sg = new SideGraph();
  _root = root;
  _referenceDupes = referenceDupes;
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
  cout << "\n-----------\n" << "ADDING " << genome->getName() << endl;
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

    set<const Genome*> inputSet;
    inputSet.insert(genome);
    inputSet.insert(candidate);
    set<const Genome*> spanningTree;
    getGenomesInSpanningTree(inputSet, spanningTree);
    if (spanningTree.size() > 1 && spanningTree.size() < dist)
    {
      dist = spanningTree.size();
      best = candidate;
    }
  }
  return best;
}

SGSequence* SGBuilder::createSGSequence(const Sequence* sequence,
                                        hal_index_t startOffset,
                                        hal_index_t length)
{
  assert(sequence != NULL);
  assert(startOffset >= 0);
  assert(length > 0);
 
  // make a new Side Graph Sequence
  stringstream name;
  name << sequence->getFullName();
  if (length < sequence->getSequenceLength())
  {
    name << "_" << startOffset << ":" << length;
  }  
  SGSequence* sgSeq = new SGSequence(-1, length, name.str());

  // add to the Side Graph
  _sg->addSequence(sgSeq);
  // keep record of where it came from (ie to trace back from the side graph
  // to hal
  _seqMapBack.insert(pair<SGSequence*, pair<const Sequence*, hal_index_t> >(
                       sgSeq,
                       pair<const Sequence*, hal_index_t>(sequence,
                                                          startOffset)));
  
  // keep record in other direction (ie to map from hal into the side graph)
  sg_seqid_t halSequenceID = (sg_seqid_t)sequence->getArrayIndex();
  assert(halSequenceID >= 0);
  cout << "ai0 ";
  _lookup->addInterval(SGPosition(halSequenceID, startOffset),
                      SGPosition(sgSeq->getID(), 0),
                      length, false);

  cout << "Create Sequence " << *sgSeq << " from hal(" << halSequenceID
       << ", " << startOffset << " l=" << length << endl;
  return sgSeq;
}

const SGJoin* SGBuilder::createSGJoin(sg_seqid_t seqId1, sg_int_t pos1,
                                      bool forward1,
                                      sg_seqid_t seqId2, sg_int_t pos2,
                                      bool forward2)
{
  SGJoin* join = new SGJoin(SGSide(SGPosition(seqId1, pos1), forward1),
                            SGSide(SGPosition(seqId2, pos2), forward2));
  return _sg->addJoin(join);
}

void SGBuilder::mapSequence(const Sequence* sequence,
                            hal_index_t globalStart,
                            hal_index_t globalEnd,
                            const Genome* target)
{
  const Genome* genome = sequence->getGenome();
  
  // first genome added: it's the reference so we add sequences
  // directly. TODO: handle self-dupes
  if (target == NULL || genome == _root || globalEnd == globalStart)
  {
    createSGSequence(sequence, globalStart - sequence->getStartPosition(),
                     globalEnd - globalStart + 1);
  }
  // only conintue if we're a second genome, or there's a possibility
  // of self dupes in the root.
  // if (genome != _root && globalEnd > globalStart)
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

    ///// DEBUG
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
    /////
    
    while (refSeg->getArrayIndex() < lastIndex &&
           refSeg->getStartPosition() <= globalEnd)  
    {
      refSeg->getMappedSegments(mappedSegments, target, &_mapPath,
                                true, 0, _root, _mapMrca);
      refSeg->toRight(globalEnd);
    }

    vector<Block*> blocks;
    blocks.reserve(mappedSegments.size());
    SGPosition prevHook(SideGraph::NullPos);

    vector<MappedSegmentConstPtr> fragments;
    BlockMapper::MSSet emptySet;
    set<hal_index_t> queryCutSet;
    set<hal_index_t> targetCutSet;
    for (set<MappedSegmentConstPtr>::iterator i = mappedSegments.begin();
         i != mappedSegments.end(); ++i)
    {
      BlockMapper::extractSegment(i, emptySet, fragments, &mappedSegments, 
                                  targetCutSet, queryCutSet);
      Block* block = new Block();
      fragmentsToBlock(fragments, *block);
      blocks.push_back(block);
    }
    // extractSegment() method above works in sorted order by target.
    // we need sorted order by source (to detect insertions, for example).
    sort(blocks.begin(), blocks.end(), BlockPtrLess());

    if (blocks.empty() == true)
    {
      // case with zero mapped blocks.  entire segment will be insertion. 
      updateSegment(NULL, NULL, NULL, prevHook, sequence,
                    genome, globalStart, globalEnd, target);
    }
    else
    {
      // we filter blocks out with cutBlock.  this offset keeps track of prev.
      size_t cutback = 0;
      for (size_t j = 0; j < blocks.size(); ++j)
      {
        Block* prev = j == 0 ? NULL : blocks[j-1-cutback];
        Block* next = j == blocks.size() - 1 ? NULL : blocks[j+1];
        Block* block = cutBlock(prev, blocks[j]);
        cutback = block == NULL ? cutback + 1 : 0;
        if (block != NULL)
        {
          updateSegment(prev, block, next, prevHook, sequence,
                        genome, globalStart, globalEnd, target);
        }
      }
    }
    for (size_t j = 0; j < blocks.size(); ++j)
    {
      delete blocks[j];
    }
  }
}


void SGBuilder::updateSegment(Block* prevBlock,
                              Block* block,
                              Block* nextBlock,
                              SGPosition& prevHook,
                              const Sequence* srcSequence,
                              const Genome* srcGenome,
                              hal_index_t globalStart, hal_index_t globalEnd,
                              const Genome* tgtGenome)
{
  // handle insertion (source sequence not mapped to target)
  hal_index_t prevSrcPos = prevBlock != NULL ? prevBlock->_srcEnd : 0;
  hal_index_t srcPos = block != NULL ? block->_srcStart :
     globalEnd - srcSequence->getEndPosition();


  cout << "PREV " << prevBlock << endl;
  cout << "CUR  " << block << endl;
  
  if (srcPos > prevSrcPos + 1)
  {
    // insert new sequence for gap between prevBock and block
    SGSequence* insertSeq = createSGSequence(srcSequence, prevSrcPos + 1,
                                              srcPos - prevSrcPos - 1);
    // join it on end of last inserted side graph position
    if (prevHook != SideGraph::NullPos)
    {
      createSGJoin(prevHook.getSeqID(), prevHook.getPos(),
                   !prevBlock->_reversed, insertSeq->getID(), 0, false);
    }

    // our new hook is the end of this new sequence
    prevHook.setSeqID(insertSeq->getID());
    prevHook.setPos(insertSeq->getLength() - 1);
  }
  if (block != NULL)
  {
    // Look up where our alignment block target (all HAL coordinates)
    // fits into the side graph. 
    SGPosition halPosition;
    halPosition.setSeqID((sg_seqid_t)block->_tgtSeq->getArrayIndex());
    halPosition.setPos(block->_reversed ? block->_tgtEnd : block->_tgtStart);
    cout << "looking up " << halPosition << endl;
    GenomeLUMap::iterator lui = _luMap.find(
      block->_tgtSeq->getGenome()->getName());
    assert(lui != _luMap.end());
    SGPosition sgPosition = lui->second->mapPosition(halPosition);
    assert(sgPosition != SideGraph::NullPos);

    // now we map in our alignment block source into the lookup table
    SGPosition newHalPosition;
    newHalPosition.setSeqID((sg_seqid_t)block->_srcSeq->getArrayIndex());
    newHalPosition.setPos(block->_srcStart);
    sg_int_t intervalLen =  block->_srcEnd - block->_srcStart + 1;
    cout << "ai1 ";
    _lookup->addInterval(newHalPosition, sgPosition, intervalLen,
                         block->_reversed);

    if (prevHook != SideGraph::NullPos)
    {
      // we can have a case where a block aligns identical segments
      // (only when looking for dupes in reference).  do not want to
      // add such self joins here
      if (isSelfBlock(*block) == false)
      {
        // we now know enough to join to prevBlock
        createSGJoin(prevHook.getSeqID(), prevHook.getPos(),
                     prevBlock ? !prevBlock->_reversed : true,
                     sgPosition.getSeqID(),
                     sgPosition.getPos(), !block->_reversed);
      }
      
    // our new hook is the end of the last join
    prevHook = sgPosition;
    }
  }
}


void SGBuilder::fragmentsToBlock(const vector<MappedSegmentConstPtr>& fragments,
                                 Block& block) const
{
  block._tgtStart = min(min(fragments.front()->getStartPosition(), 
                            fragments.front()->getEndPosition()),
                        min(fragments.back()->getStartPosition(),
                            fragments.back()->getEndPosition()));
  block._tgtEnd = max(max(fragments.front()->getStartPosition(), 
                          fragments.front()->getEndPosition()),
                      max(fragments.back()->getStartPosition(),
                          fragments.back()->getEndPosition()));
  block._tgtSeq = fragments.front()->getSequence();
  assert(block._tgtStart <= block._tgtEnd);

  SlicedSegmentConstPtr srcFront = fragments.front()->getSource();
  SlicedSegmentConstPtr srcBack = fragments.back()->getSource();

  block._srcStart =  min(min(srcFront->getStartPosition(), 
                             srcFront->getEndPosition()),
                         min(srcBack->getStartPosition(),
                             srcBack->getEndPosition()));
  block._srcEnd = max(max(srcFront->getStartPosition(), 
                          srcFront->getEndPosition()),
                      max(srcBack->getStartPosition(),
                          srcBack->getEndPosition()));
  block._srcSeq = srcFront->getSequence();
  assert(block._srcStart <= block._srcEnd);

  // convert to segment level;
  block._tgtStart -= block._tgtSeq->getStartPosition();
  block._tgtEnd -= block._tgtSeq->getStartPosition();
  block._srcStart -= block._srcSeq->getStartPosition();
  block._srcEnd -= block._srcSeq->getStartPosition();

  block._reversed = srcFront->getReversed() != fragments.front()->getReversed();
}

SGBuilder::Block* SGBuilder::cutBlock(Block* prev, Block* cur)
{
  if (prev != NULL)
  {
    //cout << "cPrev " << (size_t)prev << " " << prev << endl
    //     << "cCur " << (size_t)cur << " " << cur << endl;
    assert(prev->_srcStart <= prev->_srcEnd);
    assert(cur->_srcStart <= cur->_srcEnd);
    cur->_srcStart = max(cur->_srcStart, prev->_srcEnd + 1);
    if (cur->_srcStart > cur->_srcEnd)
    {
      cur = NULL;
    }
  }
  if (cur != NULL)
  {
    assert(cur->_srcStart <= cur->_srcEnd);
    assert(cur->_tgtStart <= cur->_tgtEnd);
  }
  //cout << "cCuO " << (size_t)cur << " " << cur << endl;

  return cur;
}
