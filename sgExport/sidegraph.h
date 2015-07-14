/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SIDEGRAPH_H
#define _SIDEGRAPH_H

/*
 * This is a quick and simple in memory side graph, with efforts
 * to use terminolgy from the latest AVRO schema doc.
 * Written expressly for hal2sg.   
 *
 */

#include <cstdlib>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

#include "sgjoin.h"
#include "sgsegment.h"

/**
 * Our simple in-memory sidegraph
 */
class SideGraph
{
public:
   // todo: replace with hash
   typedef std::set<const SGJoin*, SGJoinPtrLess> JoinSet;

   typedef std::vector<SGSequence*> SequenceSet;

   static const SGPosition NullPos;

public:
   SideGraph();
   virtual ~SideGraph();

   /**
    * Check if join is in the graph and return it else NULL.
    */
   const SGJoin* getJoin(const SGJoin* join) const;

   /**
    * Join added to the graph.  Note that SideGraph will take over
    * responsability for freeing memory of the join.
    * Duplicate joins are silently filtered and erased.  
    */
   const SGJoin* addJoin(SGJoin* join);

   /** 
    * Get the join set, ie to iterate over or something
    */
   const JoinSet* getJoinSet() const;

   /**
    * Fetch sequence using its ID
    */   
   const SGSequence* getSequence(sg_seqid_t id) const;

   /**
    * Get number of sequences
    */
   sg_int_t getNumSequences() const;
   
   /** 
    * Add a sequence into the graph.  SideGraph takes deletion resp 
    * ID field will be ignored and updated automatically
    */
   const SGSequence* addSequence(SGSequence* seq);
   
protected:

   JoinSet _joinSet;
   SequenceSet _seqSet;

private:
   SideGraph(const SideGraph& sg);
};

/** 
 * Print the Side Graph
 */
std::ostream& operator<<(std::ostream& os, const SideGraph& sg);


inline const SGJoin* SideGraph::getJoin(const SGJoin* join) const
{
  JoinSet::const_iterator i;
  // joins undirected, so we always ensure small side first
  if (join->getSide2() < join->getSide1())
  {
    SGJoin tempJoin(join->getSide2(), join->getSide1());
    i = _joinSet.find(&tempJoin);
  }
  else
  {
    i = _joinSet.find(join);
  }
  return i != _joinSet.end() ? *i : NULL;
}

inline const SGJoin* SideGraph::addJoin(SGJoin* join)
{
  assert(join->getSide1().getBase().getSeqID() >= 0);
  assert(join->getSide2().getBase().getSeqID() >= 0);
  assert(getSequence(join->getSide1().getBase().getSeqID()) != NULL);
  assert(getSequence(join->getSide2().getBase().getSeqID()) != NULL);
  // joins undirected, so we always ensure small side first
  if (join->getSide2() < join->getSide1())
  {
    join->swap();
  }
  std::pair<JoinSet::iterator, bool> r = _joinSet.insert(join);
  if (r.second == false)
  {
    delete join;
  }
  return *r.first;
}

inline const SideGraph::JoinSet* SideGraph::getJoinSet() const
{
  return &_joinSet;
}

inline const SGSequence* SideGraph::getSequence(sg_seqid_t id)
  const
{
  return _seqSet.at(id);
}

inline sg_int_t SideGraph::getNumSequences() const
{
  return _seqSet.size();
}

inline const SGSequence* SideGraph::addSequence(SGSequence* seq)
{
  seq->setID(_seqSet.size());
  _seqSet.push_back(seq);
  return seq;
}

#endif
