/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sidegraph.h"

using namespace std;

const SGPosition SideGraph::NullPos = {-1, -1};

SideGraph::SideGraph()
{

}

SideGraph::~SideGraph()
{
  for (JoinSet::iterator i = _joinSet.begin(); i != _joinSet.end(); ++i)
  {
    delete *i;
  }

  for (SequenceSet::iterator i = _seqSet.begin(); i != _seqSet.end(); ++i)
  {
    delete *i;
  }
}

