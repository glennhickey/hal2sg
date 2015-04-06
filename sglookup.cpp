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

void SGLookup::init(size_t numSequences)
{
  _mapVec = PosMapVec(numSequences);
}

void SGLookup::addInterveral(const SGPosition& inPos,
                             const SGPosition& outPos,
                             sg_int_t length,
                             bool reversed)
{
  PosMap& pm = _mapVec.at(inPos._seqid);

  sg_int_t left = inPos._pos;
  sg_int_t right = inPos._pos + length;

  // find the left point
  PosMap::iterator li = pm.insert(pair<sg_int_t, SGPosition>(
                                    left, SideGraph::NullPos)).first;

  // find the right point
  PosMap::iterator ri = li;
  ++ri;
  if (ri == pm.end() || ri->second._pos != right)
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
