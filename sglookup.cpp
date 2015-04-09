/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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
  for (size_t i = 0; i < _seqNameToId.size(); ++i)
  {
    _seqNameToId.insert(pair<string, sg_int_t>(_seqIdToName[i], sg_int_t(i)));
  }
}

void SGLookup::addInterveral(const SGPosition& inPos,
                             const SGPosition& outPos,
                             sg_int_t length,
                             bool reversed)
{
  PosMap& pm = _mapVec.at(inPos.getSeqID());

  sg_int_t left = inPos.getPos();
  sg_int_t right = inPos.getPos() + length;

  // find the left point
  PosMap::iterator li = pm.insert(pair<sg_int_t, SGPosition>(
                                    left, SideGraph::NullPos)).first;

  // find the right point
  PosMap::iterator ri = li;
  ++ri;
  if (ri == pm.end() || ri->second.getPos() != right)
  {
    ri = pm.insert(li, pair<sg_int_t, SGPosition>(
                                      right, SideGraph::NullPos));
    --ri;
    assert(ri == li);
  }

  // update the left point
  assert(li->second == SideGraph::NullPos);
  li->second = outPos;
}
