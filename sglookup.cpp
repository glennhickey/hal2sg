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

void SGLookup::addInterval(const SGPosition& inPos,
                           const SGPosition& outPos,
                           sg_int_t length,
                           bool reversed)
{
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

SGSide SGLookup::mapPosition(const SGPosition& inPos) const
{
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
  if (outSide.getForward() == false)
  {
    PosMap::const_iterator j = i;
    ++j;
    sg_int_t transform = j->first - offset - i->first - 1 - offset;
    outPos.setPos(outPos.getPos() + transform);
  }
  outSide.setBase(outPos);
  return outSide;
}
