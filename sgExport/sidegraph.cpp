/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sidegraph.h"

using namespace std;

const SGPosition SideGraph::NullPos(-1, -1);

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

ostream& operator<<(ostream& os, const SideGraph& sg)
{
  os << "SideGraph {\n";
  sg_int_t i = 0;
  for (; i < sg.getNumSequences(); ++i)
  {
    os << "Sequence " << i << ": " << *sg.getSequence(i) << "\n";
  }
  const SideGraph::JoinSet* js = sg.getJoinSet();
  i = 0;
  for (SideGraph::JoinSet::const_iterator j = js->begin(); j != js->end(); ++j)
  {
    os << "Join " << i++ << ": " << **j << "\n";
  }
  return os << "}" << endl;
}
