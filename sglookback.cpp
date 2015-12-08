/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <limits>
#include <sstream>

#include "hal.h"
#include "sglookback.h"

using namespace std;
using namespace hal;

SGLookBack::SGLookBack()
{

}

SGLookBack::~SGLookBack()
{
  
}

void SGLookBack::clear()
{
  _lookup = SGLookup();
  _halSeqMap.clear();
}

void SGLookBack::addInterval(const SGPosition& inPos, const Sequence* outSeq,
                             hal_index_t outOffset, 
                             sg_int_t length, bool reversed)
{
  if (_halSeqMap.size() <= inPos.getSeqID())
  {
    size_t oldSize = _halSeqMap.size();
    _halSeqMap.resize(inPos.getSeqID() + 1, NULL);
    size_t newSize = _halSeqMap.size();
    // little hack to resize lookup with fake names
    for (size_t i = oldSize; i < newSize; ++i)
    {
      stringstream ss;
      ss << i;
      _lookup.addSequence(ss.str());
    }
  }
  assert(_halSeqMap[inPos.getSeqID()] == outSeq ||
         _halSeqMap[inPos.getSeqID()] == NULL);
  if (_halSeqMap[inPos.getSeqID()] == NULL)
  {
    _halSeqMap[inPos.getSeqID()] = outSeq;
  }
  _lookup.addInterval(inPos, SGPosition(inPos.getSeqID(), outOffset),
                      length, reversed);
}

pair<const Sequence*, SGSide>
SGLookBack::mapPosition(const SGPosition& inPos) const
{
  assert(_halSeqMap.size() > inPos.getSeqID());
  SGSide mapSide = _lookup.mapPosition(inPos);
  return pair<const Sequence*, SGSide>(_halSeqMap[inPos.getSeqID()], mapSide);
}

void SGLookBack::getPath(const SGPosition& startPos,
                         int length,
                         bool forward,
                         vector<SGSegment>& outPath,
                         vector<const Sequence*>& outHalSeqs) const
{
  _lookup.getPath(startPos, length, forward, outPath);
  outHalSeqs.resize(outPath.size(), NULL);
  for (size_t i = 0; i < outPath.size(); ++i)
  {
    outHalSeqs[i] = _halSeqMap[outPath[i].getSide().getBase().getSeqID()];
  }
}
