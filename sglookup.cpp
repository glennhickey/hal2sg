/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <limits>

#include "sglookup.h"

using namespace std;

SGLookup::SGLookup()
{

}

SGLookup::~SGLookup()
{
  
}

void SGLookup::init(const vector<string>& sequenceNames)
{
  _mapVec = PosMapVec(sequenceNames.size());

  _seqIdToName = sequenceNames;
  _seqNameToId.clear();
  for (size_t i = 0; i < _seqIdToName.size(); ++i)
  {
    _seqNameToId.insert(pair<string, sg_int_t>(_seqIdToName[i], sg_int_t(i)));

    // add sentinel markers
    PosMap& pm = _mapVec[i];
    pm.insert(pair<sg_int_t, SGSide>(0, SGSide(SideGraph::NullPos, true)));
    pm.insert(pair<sg_int_t, SGSide>(numeric_limits<sg_int_t>::max(),
                                     SGSide(SideGraph::NullPos, true)));
  }
}

void SGLookup::addSequence(const string& sequenceName)
{
  _mapVec.push_back(PosMap());
  _seqIdToName.push_back(sequenceName);
  size_t i = _seqIdToName.size() - 1;

  _seqNameToId.insert(pair<string, sg_int_t>(_seqIdToName[i], sg_int_t(i)));

  // add sentinel markers
  PosMap& pm = _mapVec[i];
  pm.insert(pair<sg_int_t, SGSide>(0, SGSide(SideGraph::NullPos, true)));
  pm.insert(pair<sg_int_t, SGSide>(numeric_limits<sg_int_t>::max(),
                                   SGSide(SideGraph::NullPos, true)));

}

void SGLookup::addInterval(const SGPosition& inPos,
                           const SGPosition& outPos,
                           sg_int_t length,
                           bool reversed)
{
//  cout << "Add interval " << inPos << "->" << outPos << ", "
//       << length << " " << reversed << endl;
  assert(inPos.getPos() >= 0 && outPos.getPos() >= 0);
  PosMap& pm = _mapVec.at(inPos.getSeqID());

  sg_int_t left = inPos.getPos();
  sg_int_t right = inPos.getPos() + length;

  // find the left point
  SGSide leftSide(SideGraph::NullPos, true);
  PosMap::iterator li = pm.insert(pair<sg_int_t, SGSide>(
                                    left, leftSide)).first;

  // find the right point
  PosMap::iterator ri = li;
  ++ri;
  if (ri == pm.end() || ri->second.getBase().getPos() != right)
  {
    SGSide rightSide(SideGraph::NullPos, true);
    ri = pm.insert(li, pair<sg_int_t, SGSide>(right, rightSide));
    --ri;
    assert(ri == li);
  }
  
  // update the left point
  assert(li->second.getBase() == SideGraph::NullPos);
  li->second = SGSide(outPos, !reversed);
}

SGSide SGLookup::mapPosition(const SGPosition& inPos, sg_int_t* outDist,
                             bool outDistReversed) const
{
  if (outDist != NULL)
  {
    *outDist = -1;
  }
  const PosMap& pm = _mapVec.at(inPos.getSeqID());
  PosMap::const_iterator i = pm.lower_bound(inPos.getPos());
  assert(i != pm.end());

  if (i->first > inPos.getPos())
  {
    assert(i != pm.begin());
    --i;
  }
  
  if (i->second.getBase() == SideGraph::NullPos)
  {
    return i->second;
  }
  
  assert(i->first <= inPos.getPos());  
  sg_int_t offset = inPos.getPos() - i->first;
  SGSide outSide(i->second);
  SGPosition outPos(outSide.getBase());
  outPos.setPos(outPos.getPos() + offset);
  PosMap::const_iterator j = i;
  ++j;
  if (outSide.getForward() == false)
  {
    sg_int_t transform = j->first - offset - i->first - 1 - offset;
    outPos.setPos(outPos.getPos() + transform);
  }
  outSide.setBase(outPos);
  if (outDist != NULL)
  {
    // forward
    if (!outDistReversed)
    {
      *outDist = j->first - 1 - inPos.getPos();
    }
    // reversed
    else
    {
      *outDist = inPos.getPos() - i->first;
    }
    assert(outDist >= 0);
  }
  return outSide;
}

void SGLookup::getPath(const SGPosition& startPos,
                       const SGPosition& endPos,
                       vector<SGSegment>& outPath) const
{
  SGPosition halStart = startPos;
  SGPosition halEnd = endPos;
  bool backward = endPos < startPos;
  if (backward == true)
  {
    // always make a forward path.  if query is reversed, we will
    // remember here and flip at very end. 
    swap(halStart, halEnd);
  }
  
  assert(halStart.getSeqID() == halEnd.getSeqID());

  const PosMap& pm = _mapVec.at(halStart.getSeqID());

  // find marker to the right of halStart
  PosMap::const_iterator i = pm.lower_bound(halStart.getPos());
  assert(i != pm.end());
  if (i->first == halStart.getPos())
  {
    ++i;
  }
  // find marker >= halEnd (ie one-past what we iterate)
  PosMap::const_iterator j = pm.lower_bound(halEnd.getPos());
  //assert(j != pm.begin());
  if (j->first == halEnd.getPos())
  {
    ++j;
  }

  outPath.clear();
  sg_int_t pathLength = 0;

  sg_int_t prevHalPos = halStart.getPos();
  SGSide prevSgSide = mapPosition(halStart);

  for (PosMap::const_iterator k = i; k != j; ++k)
  {
    assert(k != pm.end());
    assert(k != pm.begin());
    sg_int_t halPos = k->first;
    sg_int_t segLen = halPos - prevHalPos;
    assert(segLen > 0);
    SGSide sgSide = k->second;

    // note we are taking advantage of the fact that
    // the side returned by mapPosition has a forward flag
    // consistent with its use in sgsegment.
    if (prevSgSide.getForward() == false && k != i)
    {
      // interval lookup always stores intervals left->right.
      // for reverse mapping, we have to manually adjust segment
      // coordinate (unless first iteration, which is already adjusted
      // by call to mapPosition)
      prevSgSide.setBase(SGPosition(prevSgSide.getBase().getSeqID(),
                                    prevSgSide.getBase().getPos() + segLen -1));
    }
    SGSegment nextSeg(prevSgSide, segLen);
    // merge two consecutive segments that can be walked without a join
    if (!outPath.empty() &&
        outPath.back().getSide().getForward() ==
        nextSeg.getSide().getForward() &&
        outPath.back().getOutSide().lengthTo(nextSeg.getInSide()) == 0)
    {
      outPath.back().setLength(outPath.back().getLength() +
                               nextSeg.getLength());
    }
    else
    {
      outPath.push_back(SGSegment(prevSgSide, segLen));
    }
    pathLength += segLen;
    prevHalPos = halPos;
    prevSgSide = sgSide;
    assert(outPath.back().getMinPos().getPos() >=0);
  }

  sg_int_t segLen = halEnd.getPos() - prevHalPos + 1;
  assert(segLen > 0);
  if (prevSgSide.getForward() == false
      // gack, only need to flip if loop entered.  this is ugly but fixes:
      && i != j)
  {
    prevSgSide.setBase(SGPosition(prevSgSide.getBase().getSeqID(),
                                  prevSgSide.getBase().getPos() + segLen -1));
  }
  outPath.push_back(SGSegment(prevSgSide, segLen));
  pathLength += segLen;
  (void)pathLength;
  assert(pathLength == halEnd.getPos() - halStart.getPos() + 1);

  // we really wanted our path in the other direction.  flip the
  // order of the vector, and the orientation of every segment. 
  if (backward == true)
  {
    reverse(outPath.begin(), outPath.end());
    for (size_t i = 0; i < outPath.size(); ++i)
    {
      outPath[i].flip();
    }
  }  
}
