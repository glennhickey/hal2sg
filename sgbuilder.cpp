/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <sstream>

#include "sgbuilder.h"
#include "snphandler.h"
#include "halBlockMapper.h"


using namespace std;
using namespace hal;

SGBuilder::SGBuilder() : _sg(0), _root(0), _mapRoot(0), _lookup(0), _mapMrca(0),
                         _referenceDupes(true), _camelMode(false),
                         _snpHandler(0),
                         _onlySequenceNames(false)
{

}

SGBuilder::~SGBuilder()
{
  clear();
}

void SGBuilder::init(AlignmentConstPtr alignment, const Genome* root,
                     bool referenceDupes, bool camelMode,
                     bool onlySequenceNames)
{
  clear();
  _alignment = alignment;
  _sg = new SideGraph();
  _root = root;
  _mapRoot = root;
  _referenceDupes = referenceDupes;
  _camelMode = camelMode;
  _snpHandler = new SNPHandler(_sg, false, onlySequenceNames);
  _onlySequenceNames = onlySequenceNames;
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
  _lookBack.clear();
  _mapPath.clear();
  _mapMrca = NULL;
  _firstGenomeName.erase();
  _halSequences.clear();
  delete _snpHandler;
  _snpHandler = NULL;
}

const SideGraph* SGBuilder::getSideGraph() const
{
  return _sg;
}

size_t SGBuilder::getSequenceString(const SGSequence* sgSequence,
                                    string& outString,
                                    sg_int_t pos,
                                    sg_int_t length) const
{
  outString.clear();
  hal_index_t len = length == -1 ? sgSequence->getLength() : length;
  vector<SGSegment> segPath;
  vector<const Sequence*> halSeqPath;
  _lookBack.getPath(SGPosition(sgSequence->getID(), pos),
                    SGPosition(sgSequence->getID(), pos + len - 1),
                    segPath, halSeqPath);
  for (size_t i = 0; i < segPath.size(); ++i)
  {
    const Sequence* halSeq = halSeqPath[i];
    const SGSegment& seg = segPath[i];
    sg_int_t leftCoord = seg.getMinPos().getPos();
    string buffer;
    if (_camelMode == true && halSeq->getGenome()->getParent() == NULL)
    {
      // adams root has not sequence. we infer it from children as a hack.
      getRootSubString(buffer, halSeq, leftCoord, seg.getLength());
    }
    else
    {
      halSeq->getSubString(buffer, leftCoord, seg.getLength());
    }
    if (seg.getSide().getForward() == false)
    {
      reverseComplement(buffer);
    }
    outString += buffer;
  }
  assert(outString.length() == len);
  return outString.length();
}

const string& SGBuilder::getPrimaryGenomeName() const
{
  return _firstGenomeName;
}

string SGBuilder::getHalGenomeName(const SGSequence* sgSequence) const
{
  return _lookBack.getHalGenomeName(sgSequence);
}

const vector<const Sequence*>& SGBuilder::getHalSequences() const
{
  return _halSequences;
}

void SGBuilder::getHalSequencePath(const Sequence* halSeq,
                                   vector<SGSegment>& outPath) const
{
  GenomeLUMap::const_iterator lui = _luMap.find(halSeq->getGenome()->getName());
  assert(lui != _luMap.end());
  SGPosition start(halSeq->getArrayIndex(), 0);
  SGPosition end(halSeq->getArrayIndex(),
                 max((hal_size_t)0, halSeq->getSequenceLength() - 1));
  lui->second->getPath(start, end, outPath);
}

void SGBuilder::addGenome(const Genome* genome,
                          const Sequence* sequence, 
                          hal_index_t start,
                          hal_index_t length)
{
  assert(_luMap.find(genome->getName()) == _luMap.end());
  if (_firstGenomeName.empty())
  {
    // make note of reference genome to tell if sequences are
    // derived or not in sql. 
    _firstGenomeName = genome->getName();
  }

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
    for (; si != sie; si->toNext())
    {
      const Sequence* curSequence = si->getSequence();
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
    // only path though root when self aligning. 
    const Genome* pathStart = genome == target ? _root : genome;
    inputSet.insert(pathStart);
    inputSet.insert(target);
    getGenomesInSpanningTree(inputSet, _mapPath);
    _mapRoot = getLowestCommonAncestor(_mapPath);
  }

  // Convert sequence by sequence
  for (size_t i = 0; i < seqNames.size(); ++i)
  {
    const Sequence* curSequence = genome->getSequence(seqNames[i]);
    hal_index_t endPos = start + length - 1;
    hal_index_t curStart = std::max(start, curSequence->getStartPosition());
    hal_index_t curEnd = std::min(endPos, curSequence->getEndPosition());

    // note all coordinates global
    if (curStart <= curEnd)
    {
      mapSequence(curSequence, curStart, curEnd, target);
      _halSequences.push_back(curSequence);
    }
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

pair<SGSide, SGSide> SGBuilder::createSGSequence(const Sequence* sequence,
                                                 hal_index_t startOffset,
                                                 hal_index_t length)
{
  assert(sequence != NULL);
  assert(startOffset >= 0);
  assert(length > 0);

  // compute Dupes
  vector<Block*> rawBlocks;
  if (sequence->getGenome() != _mapRoot)
  {
    cout << "compute blocks " << sequence->getName() << " "
         << startOffset  << ": " << length << endl;
    computeBlocks(sequence, sequence->getStartPosition() + startOffset,
                  sequence->getStartPosition() + startOffset + length - 1,
                  NULL, rawBlocks);
  }

  // clean up weird overlap cases from self-dupe code (to revisit)
  vector<Block*> blocks;
  filterOverlaps(rawBlocks, blocks);
  if (blocks.size() > 0)
  {
    cout << "Dupe blocks found = ";
    for (size_t i = 0; i < blocks.size(); ++i)
    {
      cout << blocks[i] << ", ";
    }
    cout << endl;
  }

  // for each block, a flag if it's collapsed or not
  vector<bool> collapsed(blocks.size(), false);
  getCollapsedFlags(blocks, collapsed);  

  // We can optimize this out probably but one pass just to compute
  // the length of the new sequence
  hal_index_t newSeqLen = 0;
  hal_index_t prev = -1;
  for (size_t i = 0; i < blocks.size(); ++i)
  {
    Block* block = blocks[i];
    hal_index_t delta = block->_srcStart - prev;
    // add gap between blocks
    if (delta > 1)
    {
      newSeqLen += delta - 1;
    }
    if (collapsed[i] == false)
    {
      newSeqLen += block->_srcEnd - block->_srcStart + 1;
    }
    prev = block->_srcEnd;
  }
  if (length -1 > prev)
  {
    newSeqLen += length - prev - 1;
  }
  assert(newSeqLen <= length);

  // make a new Side Graph Sequence
  stringstream name;
  name << getHalSeqName(sequence);
  if (length < (hal_index_t)sequence->getSequenceLength())
  {
    name << "_" << startOffset << "_" << newSeqLen;
  }  
  SGSequence* sgSeq = new SGSequence(-1, newSeqLen, name.str());

  // add to the Side Graph
  _sg->addSequence(sgSeq);

  // update the _lookup structure for the entire new sequence
  // (gaps and uncollapsed regions)
  updateDupeBlockLookups(sequence, startOffset, length, blocks, collapsed,
                         sgSeq);
  
  // we have all pairwise alignment blocks in our list.  we only
  // need blocks where the target range is uncollapsed.  filter
  // everything else out here.
  filterRedundantDupeBlocks(blocks, sgSeq);

  // add all the joins and snps resulting from the duplication blocks
  addDupeJoins(sequence, startOffset, length, blocks, sgSeq);
    
  // Hooks to the start and end of the input subsequence in the sidegraph.
  // note to self: these two maps could conceivably be optimized out.
  pair<SGSide, SGSide> outHooks;
  sg_seqid_t halSequenceID = (sg_seqid_t)sequence->getArrayIndex();
  outHooks.first =  _lookup->mapPosition(SGPosition(halSequenceID,
                                                    startOffset));
  outHooks.second = _lookup->mapPosition(SGPosition(halSequenceID,
                                                    startOffset + length-1));
  // if not reversed : false
  outHooks.first.setForward(!outHooks.first.getForward());
  // if not revered : true
  outHooks.second.setForward(outHooks.second.getForward());
  assert(outHooks.first.getBase() != SideGraph::NullPos);
  assert(outHooks.second.getBase() != SideGraph::NullPos);

  return outHooks;
}

const SGJoin* SGBuilder::createSGJoin(const SGSide& side1, const SGSide& side2)
{
  SGJoin* join = new SGJoin(side1, side2);

  // filter trivial join
/*  if (abs(join->getSide1().getBase().getPos() -
    join->getSide2().getBase().getPos()) == 1 &&
    join->getSide1().getForward() !=
    join->getSide2().getForward())
    {
    cout << "filter trivial " << *join << endl;
    delete join;
    return NULL;
    }
*/
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
  if (target == NULL || genome == _root)
  {
    createSGSequence(
      sequence, globalStart - sequence->getStartPosition(),
      globalEnd - globalStart + 1);
  }
  // only conintue if we're a second genome, or there's a possibility
  // of self dupes in the root.
  // if (genome != _root && globalEnd > globalStart)
  else 
  {
    vector<Block*> blocks;
    computeBlocks(sequence, globalStart, globalEnd, target, blocks);

    SGSide prevHook(SideGraph::NullPos, true);
    _lastJoin = NULL;
   
    // we will deal in sequence-relative (as opposed to hals wonky
    // genome-relative) coordinates from here on.
    hal_index_t sequenceStart = globalStart - sequence->getStartPosition();
    hal_index_t sequenceEnd = globalEnd - sequence->getStartPosition();
    
    if (blocks.empty() == true)
    {
      // case with zero mapped blocks.  entire segment will be insertion. 
      visitBlock(NULL, NULL, NULL, prevHook, sequence,
                   genome, sequenceStart, sequenceEnd, target);
    }
    else
    {
      // we filter blocks out with cutBlock.  this offset keeps track of prev.
      // todo: we need to put BlockMapper machinery to cut and self-map
      // to handle the reference/insertion paralogy case (maybe).
      size_t cutback = 0;
      for (size_t j = 0; j < blocks.size(); ++j)
      {
        Block* prev = j == 0 ? NULL : blocks[j-1-cutback];
        Block* next = j == blocks.size() - 1 ? NULL : blocks[j+1];
        Block* block = cutBlock(prev, blocks[j]);

        cutback = block == NULL ? cutback + 1 : 0;
        if (block != NULL)
        {
          visitBlock(prev, block, next, prevHook, sequence,
                       genome, sequenceStart, sequenceEnd, target);
        }
      }
      // add insert at end / last step in path
      visitBlock(blocks.back(), NULL, NULL, prevHook, sequence,
                   genome, sequenceStart, sequenceEnd, target);
    }
    for (size_t j = 0; j < blocks.size(); ++j)
    {
      delete blocks[j];
    }

    if (_lastJoin != NULL)
    {
      _sg->addJoin(_lastJoin);
      _lastJoin = NULL;
    }
  }
}

void SGBuilder::computeBlocks(const Sequence* sequence,
                              hal_index_t globalStart,
                              hal_index_t globalEnd,
                              const Genome* target,
                              vector<Block*>& blocks)
{
  const Genome* genome = sequence->getGenome();

  const Genome* currentMRCA = _mapMrca;
  if (target == NULL)
  {
    // map to self
    target = genome;
    currentMRCA = genome;
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
  cout << "Mapping " << getHalSeqName(sequence) << " to "
       << target->getName()
       << "; MRC = " << currentMRCA->getName()
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
                              true, 0, _mapRoot, currentMRCA);
    refSeg->toRight(globalEnd);
  }

  blocks.clear();
  blocks.reserve(mappedSegments.size());

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
    // filter trivial self alignments while at it
    if (block->_srcSeq == block->_tgtSeq &&
         block->_srcStart == block->_tgtStart &&
         block->_srcEnd == block->_tgtEnd)
    {
      delete block;
      block = NULL;
    }
    else
    {
      blocks.push_back(block);
    }
  }

  // extractSegment() method above works in sorted order by target.
  // we need sorted order by source (to detect insertions, for example).
  
  sort(blocks.begin(), blocks.end(), BlockPtrLess());
}

void SGBuilder::visitBlock(Block* prevBlock,
                           Block* block,
                           Block* nextBlock,
                           SGSide& prevHook,
                           const Sequence* srcSequence,
                           const Genome* srcGenome,
                           hal_index_t sequenceStart,
                           hal_index_t sequenceEnd,
                           const Genome* tgtGenome)
{
  hal_index_t prevSrcPos;  // global hal coord of end of last block
  hal_index_t srcPos;  // global hal coord of beginning of block
  if (prevBlock == NULL)
  {
    // minus one here because we consider prevSrcPos covered and want
    // to insert stuff that comes after. 
    prevSrcPos = sequenceStart - 1;
  }
  else
  {
    prevSrcPos = prevBlock->_srcEnd;
  }
  if (block != NULL)
  {
    srcPos = block->_srcStart;
  }
  else
  {
    srcPos = sequenceEnd + 1;
  }

  cout << endl << "PREV " << prevBlock << endl;
  cout << "CUR  " << block << endl;
  cout << "prevSrcPos " << prevSrcPos << " srcPos " << srcPos << endl;

  if (prevBlock == 0)
  {
    _pathLength = 0;
  }

  if (srcPos > prevSrcPos + 1)
  {    
    // handle insertion (source sequence not mapped to target)
    // insert new sequence for gap between prevBock and block
    pair<SGSide, SGSide> seqHooks = createSGSequence(srcSequence,
                                                     prevSrcPos + 1, 
                                                     srcPos - prevSrcPos - 1);
    
    // join it on end of last inserted side graph position
    if (prevHook.getBase() != SideGraph::NullPos)
    {
      createSGJoin(prevHook, seqHooks.first);
    }

    // our new hook is the end of this new sequence
    prevHook = seqHooks.second;
    _pathLength += srcPos - prevSrcPos - 1;
  }
  if (block != NULL)
  {
    mapBlockEnds(block, prevHook);
  }
}

void SGBuilder::mapBlockEnds(Block* block,
                             SGSide& prevHook)
{
  assert(block != NULL);
  
  sg_int_t blockLength = block->_srcEnd - block->_srcStart + 1;

  // We want to find the two sides in our Side Graph where the current
  // block's end points map to.  We have to be careful because it's a
  // 2-step mapping (block.src -> block.tgt -> side graph), so need to
  // account for double-reversal case.
  pair<SGSide, SGSide> blockEnds;

  // map from src to target (first mapping);
  SGPosition halTgtFirst(
    (sg_seqid_t)block->_tgtSeq->getArrayIndex(),
    block->_reversed ? block->_tgtEnd : block->_tgtStart);
    
  // map from tgt to side graph (second mapping);
  GenomeLUMap::iterator lui = _luMap.find(
    block->_tgtSeq->getGenome()->getName());
  assert(lui != _luMap.end());
  blockEnds.first = lui->second->mapPosition(halTgtFirst);
  bool sgMapReversed = !blockEnds.first.getForward();
  blockEnds.second = blockEnds.first;
  SGPosition secondPos = blockEnds.second.getBase();
  if (sgMapReversed != block->_reversed)
  {
    secondPos.setPos(secondPos.getPos() - blockLength + 1);
    blockEnds.second.setBase(secondPos);
    blockEnds.first.setForward(true);
    blockEnds.second.setForward(false);
  }
  else
  {
    secondPos.setPos(secondPos.getPos() + blockLength - 1);
    blockEnds.second.setBase(secondPos);
    blockEnds.first.setForward(false);
    blockEnds.second.setForward(true);
  }

  // we've found how our block fits into the graph (blockEnds).  Now
  // we need to map the inside of the block.  this means processing
  // all the snps, as well as making sure the lookup structure is
  // updated.
  pair<SGSide, SGSide> mappedBlockEnds =  mapBlockBody(block, blockEnds);
 
  // cout << "BE " << blockEnds.first << ", " << blockEnds.second << endl;
  // cout << "TC " << mappedBlockEnds.first << ", " << mappedBlockEnds.second
  //      << endl;

  _pathLength += blockLength;
    
  if (prevHook.getBase() != SideGraph::NullPos)
  {
    // we can have a case where a block aligns identical segments
    // (only when looking for dupes in reference).  do not want to
    // add such self joins here
    if (isSelfBlock(*block) == false)
    {
      // we now know enough to join to prevBlock
      createSGJoin(prevHook, mappedBlockEnds.first);
    }
  }
      
  // our new hook is the end of the last join
  prevHook = mappedBlockEnds.second;
      
  // note to self, this shoudl be in createSGJoin or
  // somehow centralized. 
  if (_lastJoin != NULL)
  {
    _sg->addJoin(_lastJoin);
    _lastJoin = NULL;
  }
}

pair<SGSide, SGSide>
SGBuilder::mapBlockBody(const Block* block,
                        const pair<SGSide, SGSide>& blockEnds)
{
  // note to self:  block is the pairwise HAL alignment
  //                blockEnds are the endpoints in the Side Graph
  hal_index_t length = block->_srcEnd - block->_srcStart + 1;
  string srcDNA;
  string tgtDNA;
  block->_srcSeq->getSubString(srcDNA, block->_srcStart, length);
  block->_tgtSeq->getSubString(tgtDNA, block->_tgtStart, length);
  if (block->_reversed == true)
  {
    reverseComplement(tgtDNA);
  }
  pair<SGSide, SGSide> outBlockEnds = blockEnds;

  SGSide prevHook(blockEnds.first);
  bool sgForwardMap = blockEnds.first <= blockEnds.second; 
  bool runningSnp = false;

  hal_index_t srcPos;
  hal_index_t tgtPos;
  char srcVal;
  char tgtVal;
  hal_index_t bp = 0;

  // slice block into runs of consecutive SNPs and the regions in between
  // ex:  AACGTATAC
  //      ACCGTGGAG
  // would translate to 5 slices: 123334455
  for (hal_index_t i = 0; i < length; ++i)
  {
    srcPos =  block->_srcStart + i;
    srcVal = srcDNA[i];
    tgtVal = tgtDNA[i];
    tgtPos = !block->_reversed ?  block->_tgtStart + i : block->_tgtEnd;
    bool snp = !_camelMode && _snpHandler->isSub(srcVal, tgtVal);

    if (i > 0 && snp != runningSnp)
    {
      pair<SGSide, SGSide> sliceEnds = mapBlockSlice(block, blockEnds,
                                                     bp, i - 1,
                                                     runningSnp, sgForwardMap,
                                                     srcDNA, tgtDNA);
      if (bp == 0)
      {
        outBlockEnds.first = sliceEnds.first;
      }
      else
      {
        createSGJoin(prevHook, sliceEnds.first);
      }
      prevHook = sliceEnds.second;
      bp = i;
    }
    if (i == length - 1)
    {
      pair<SGSide, SGSide> sliceEnds = mapBlockSlice(block, blockEnds, bp, i,
                                                     snp, sgForwardMap,
                                                     srcDNA, tgtDNA);
      if (bp == 0)
      {
        outBlockEnds.first = sliceEnds.first;
      }
      else
      {
        createSGJoin(prevHook, sliceEnds.first);
      }
      outBlockEnds.second = sliceEnds.second;
      prevHook = sliceEnds.second;
    }
    
    runningSnp = snp;
  }

  return outBlockEnds;
}

pair<SGSide, SGSide>
SGBuilder::mapBlockSlice(const Block* block,
                         const pair<SGSide, SGSide>& blockEnds,
                         hal_index_t srcStartOffset,
                         hal_index_t srcEndOffset,
                         bool snp, 
                         bool sgForwardMap,
                         const string& srcDNA,
                         const string& tgtDNA)
{
  SGPosition srcHalPosition((sg_seqid_t)block->_srcSeq->getArrayIndex(),
                            block->_srcStart + srcStartOffset);
  SGPosition tgtHalPosition((sg_seqid_t)block->_tgtSeq->getArrayIndex(),
                            block->_tgtStart + srcStartOffset);
  sg_int_t blockLength = srcEndOffset - srcStartOffset + 1;

  pair<SGSide, SGSide> outEnds = blockEnds;
  bool reversed = blockEnds.second < blockEnds.first;

  SGPosition startPos = blockEnds.first.getBase();
  SGPosition endPos = blockEnds.second.getBase();
  if (reversed == false)
  {
    startPos.setPos(startPos.getPos() + srcStartOffset);
    endPos.setPos(startPos.getPos() + blockLength - 1);
  }
  else
  {
    startPos.setPos(startPos.getPos() - srcStartOffset);
    endPos.setPos(startPos.getPos() - blockLength + 1);
  }
  // cout << "mapBlockSlice(" << block << "\n"
  //      << blockEnds.first << ", " << blockEnds.second << "\n"
  //      << srcStartOffset << ", " << srcEndOffset << "\n"
  //      << snp << "\n"
  //      << sgForwardMap << "\n" << endl;
  if (snp == false)
  {
    outEnds.first.setBase(startPos);
    outEnds.second.setBase(endPos);
    _lookup->addInterval(srcHalPosition, (reversed ? endPos : startPos),
                         blockLength, block->_reversed);
  }
  else
  {
    outEnds = _snpHandler->createSNP(srcDNA,
                                     tgtDNA,
                                     srcStartOffset,
                                     blockLength,
                                     block->_srcSeq,
                                     srcHalPosition,
                                     (reversed ? endPos : startPos),
                                     reversed,
                                     _lookup,
                                     &_lookBack);
  }

  return outEnds;
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

SGBuilder::Block* SGBuilder::cutBlock(Block* prev, Block* cur,
                                      bool leaveExactOverlaps)
{
  if (prev != NULL)
  {
    //cout << "cPrev " << (size_t)prev << " " << prev << endl
    //     << "cCur " << (size_t)cur << " " << cur << endl;
    assert(prev->_srcStart <= prev->_srcEnd);
    assert(cur->_srcStart <= cur->_srcEnd);
    sg_int_t cutLen = prev->_srcEnd - cur->_srcStart + 1;
    if (cutLen > 0 && (!leaveExactOverlaps ||
                       cutLen < cur->_srcEnd - cur->_srcStart + 1))
    {
      assert(false);
      cur->_srcStart += cutLen;
      cur->_tgtStart += cutLen;
      
      if (cur->_srcStart > cur->_srcEnd)
      {
        cur = NULL;
      }
    }
  }
  if (cur != NULL)
  {
    assert(cur->_srcStart <= cur->_srcEnd);
    assert(cur->_tgtStart <= cur->_tgtEnd);
  }
  return cur;
}

bool SGBuilder::verifyPath(const Sequence* sequence,
                           const vector<SGSegment>& path) const
{
  string pathString;
  string buffer;
  for (size_t i = 0; i < path.size(); ++i)
  {
    const SGSegment& seg = path[i];
    const SGSequence* seq = _sg->getSequence(
      seg.getSide().getBase().getSeqID());
    getSequenceString(seq, buffer, seg.getMinPos().getPos(), seg.getLength());
    if (seg.getSide().getForward() == false)
    {
      reverseComplement(buffer);
    }
    pathString.append(buffer);
    if (i > 0)
    {
      SGSide srcSide = path[i-1].getOutSide();
      SGSide tgtSide = seg.getInSide();
      SGJoin queryJoin(srcSide, tgtSide);      
      const SGJoin* join = _sg->getJoin(&queryJoin);
      assert(join != NULL);
      if (join == NULL)
      {
        return false;
      }
    }
  }

  if (_camelMode == true && sequence->getGenome()->getParent() == NULL)
  {
    getRootSubString(buffer, sequence, 0, sequence->getSequenceLength());
  }
  else
  {
    sequence->getString(buffer);
  }

  if (buffer != pathString)
  {
    cout << getHalSeqName(sequence) << endl;
    cout << "BUF " << buffer.length()   << endl;
    cout << "PAT " << pathString.length()  << endl;

    cout << "BUF " << buffer   << endl;
    cout << "PAT " << pathString  << endl;


    for (size_t x = 0; x < buffer.length(); ++x)
    {
      if (buffer[x] != pathString[x])
      {
        cout << x << " " << buffer[x] << "->" << pathString[x] << endl;
        break;
      }
    }
    cout << endl;
  }
  assert(buffer == pathString);
  if (buffer != pathString)
  {
    return false;
  }

  return true;
}


void SGBuilder::getRootSubString(string& outDNA, const Sequence* sequence,
                                 hal_index_t pos, hal_index_t length) const
{
  string buffer;
  outDNA.erase();
  // copy the DNA sequence up from the children for the entire genome. 
  if (_rootString.length() == 0)
  {
    const Genome* root = sequence->getGenome();
    assert (root->getParent() == NULL);
    BottomSegmentIteratorConstPtr bottom = root->getBottomSegmentIterator();
    TopSegmentIteratorConstPtr top = root->getChild(0)->getTopSegmentIterator();
    for (size_t i = 0 ; i < root->getNumBottomSegments(); ++i)
    {
      for (size_t j = 0; j < root->getNumChildren(); ++j)
      {
        if (bottom->hasChild(j))
        {
          top->toChild(bottom, j);
          top->getString(buffer);
          _rootString.append(buffer);
          break;
        }
        assert(j != root->getNumChildren() - 1);
      }
      bottom->toRight();
    }
    assert(root->getSequenceLength() == _rootString.length());
  }

  // do our query on this new root.
  outDNA = _rootString.substr(sequence->getStartPosition() + pos, length);
}

void SGBuilder::filterOverlaps(const vector<Block*>& rawBlocks,
                               vector<Block*>& blocks)
{
  // filter out overlaps.  don't think they ever happens but
  // need to revisit this logic.  Will certainly never come into
  // play on a star tree, or cases where we don't cross crazy distance
  // between output species (can restrict at higher level).  
  Block* prev = NULL;
  for (size_t i = 0; i < rawBlocks.size(); ++i)
  {
    Block* cur = cutBlock(prev, rawBlocks[i], true);
    if (cur == NULL)
    {
      delete rawBlocks[i];
    }
    else
    {
      blocks.push_back(cur);
      prev = cur;
    }
  }
}

void SGBuilder::getCollapsedFlags(const vector<Block*>& blocks,
                                  vector<bool>& collapseBlock)
{
  // a temporary interval map, where we ignore sequence id.  
  SGLookup locLook;
  locLook.init(vector<string>(1, ""));
  
  // for each block, a flag if it's collapsed or not
  for (size_t i = 0; i < blocks.size(); ++i)
  {
    Block* block = blocks[i];
    hal_index_t blockLength = block->_srcEnd - block->_srcStart + 1;
    SGPosition tgtHalPosition((sg_seqid_t)block->_tgtSeq->getArrayIndex(),
                              block->_tgtStart);
    // check if it's a duplication in a *different* side graph sequence
    SGSide mapResult = _lookup->mapPosition(tgtHalPosition);
    bool visited = mapResult.getBase() != SideGraph::NullPos;
    // check if it's a duplication in the same side graph sequence
    mapResult = locLook.mapPosition(SGPosition(0, block->_tgtStart));
    bool tgtVis = mapResult.getBase() != SideGraph::NullPos;
    mapResult = locLook.mapPosition(SGPosition(0, block->_srcStart));
    bool srcVis = mapResult.getBase() != SideGraph::NullPos;

    collapseBlock[i] = visited || tgtVis || srcVis;
    if (!tgtVis)
    {
      locLook.addInterval(SGPosition(0, block->_tgtStart),
                          SGPosition(0, block->_tgtStart),
                          blockLength, false);
    }
    if (!srcVis)
    {
      locLook.addInterval(SGPosition(0, block->_srcStart),
                          SGPosition(0, block->_srcStart),
                          blockLength, false);      
    }
  }
}

void SGBuilder::updateDupeBlockLookups(const Sequence* sequence,
                                       hal_index_t startOffset,
                                       hal_index_t length,
                                       const vector<Block*>& blocks,
                                       const vector<bool>& collapsed,
                                       const SGSequence* sgSeq)
{
  // 1 pass just to update map structures for uncollapsed
  sg_seqid_t halSequenceID = (sg_seqid_t)sequence->getArrayIndex();
  assert(halSequenceID >= 0);
  hal_index_t prev = -1;
  hal_index_t newSeqLen = 0;
  for (size_t i = 0; i < blocks.size(); ++i)
  {
    Block* block = blocks[i];
    hal_index_t delta = block->_srcStart - prev;

    if (delta > 1)
    {
      // add gap between blocks [prev+1 -> srcStart-1]
      _lookup->addInterval(SGPosition(halSequenceID, startOffset + prev + 1),
                           SGPosition(sgSeq->getID(), newSeqLen),
                           delta -1, false);  
      newSeqLen += delta - 1;
    }
    if (collapsed[i] == false)
    {
      // add map of block source to sidegraph sequence
      _lookup->addInterval(SGPosition(halSequenceID, block->_srcStart),
                           SGPosition(sgSeq->getID(), newSeqLen),
                           block->_srcEnd - block->_srcStart + 1, false);

      
      newSeqLen += block->_srcEnd - block->_srcStart + 1;
    }
    prev = block->_srcEnd;
  }
  if (length -1 > prev)
  {
    // add gap between blocks [prev+1 -> end]
    _lookup->addInterval(SGPosition(halSequenceID, startOffset + prev + 1),
                         SGPosition(sgSeq->getID(), newSeqLen),
                         length - prev -1, false);  
    newSeqLen += length - prev - 1;
  }

  // NEED TO FIX::
  // keep record of where it came from (ie to trace back from the side graph
  // to hal
  _lookBack.addInterval(SGPosition(sgSeq->getID(), 0), sequence, startOffset,
                       length, false);
//  _seqMapBack.insert(
//    pair<const SGSequence*, pair<const Sequence*, hal_index_t> >(
//      sgSeq, pair<const Sequence*, hal_index_t>(sequence, startOffset)));

}

void SGBuilder::filterRedundantDupeBlocks(vector<Block*>& blocks,
                                          const SGSequence* sgSeq)
{
  vector<Block*> filteredBlocks;
  for (size_t i = 0; i < blocks.size(); ++i)
  {
    if (_lookup->mapPosition(SGPosition(0, blocks[i]->_tgtStart)).getBase() !=
        SideGraph::NullPos)
    {
      filteredBlocks.push_back(blocks[i]);
    }
  }
  swap(filteredBlocks, blocks);
}

void SGBuilder::addDupeJoins(const Sequence* sequence,
                             hal_index_t startOffset,
                             hal_index_t length,
                             const vector<Block*>& blocks,                
                             const SGSequence* sgSeq)
{

}
