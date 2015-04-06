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

typedef size_t sg_size_t;
typedef long sg_int_t;
typedef sg_size_t sg_seqid_t;
typedef sg_size_t sg_joinid_t;

/**
 * Named sequence.  We will rely on separate map to/from HAL sequences
 * in order to get DNA, I think...  IDs will be unique numbers from
 * 0,1,2... since no plans to do anything but add, at this point.  
 */
struct SGSequence
{
   sg_seqid_t _id;
   sg_int_t _length;
   std::string _name;
};

/**
 *  Position in the sidegraph
 */
struct SGPosition
{
   sg_seqid_t _seqid;
   sg_int_t _pos;
};

inline bool operator<(const SGPosition& p1, const SGPosition& p2) {
  return p1._seqid < p2._seqid || (p1._seqid == p2._seqid &&
                                   p1._pos < p2._pos);
}

inline bool operator==(const SGPosition& p1, const SGPosition& p2) {
  return p1._seqid == p2._seqid && p1._pos == p2._pos;
}

inline bool operator!=(const SGPosition& p1, const SGPosition& p2) {
  return !(p1 == p2);
}

inline std::ostream& operator<<(std::ostream& os,
                                const SGPosition& p){
  return os << "p(" << p._seqid << "," << p._pos << ")";
}

/**
 * Position of join endpoint
 */ 
struct SGSide
{
   SGPosition _base;
   bool _forward;
};

inline bool operator<(const SGSide& s1, const SGSide& s2) {
  return s1._base < s2._base || (s1._forward == s2._forward &&
                                 s1._base < s2._base);
}

inline bool operator==(const SGSide& s1, const SGSide& s2) {
  return s1._base == s2._base && s1._forward == s2._forward;
}

inline std::ostream& operator<<(std::ostream& os, const SGSide& s) {
  return os << "s(" << s._base << "," << s._forward << ")";
}


/**
 * Join from one side to another
 */
struct Join
{
   SGSide _side1;
   SGSide _side2;
};

inline bool operator<(const Join& j1, const Join& j2) {
  return j1._side1 < j2._side1 || (j1._side1 == j2._side1 &&
                                 j1._side2 < j2._side2);
}

inline bool operator==(const Join& j1, const Join& j2) {
  return j1._side1 == j2._side1 && j1._side2 == j2._side2;
}

struct JoinPtrLess {
   bool operator()(const Join* j1, const Join* j2) const {
     return *j1 < *j2;
   }
};

inline std::ostream& operator<<(std::ostream& os, const Join& j) {
  return os << "j(" << j._side1 << "," << j._side2 << ")";
}

/**
 * Our simple in-memory sidegraph
 */
class SideGraph
{
public:
   // todo: replace with hash
   typedef std::set<const Join*, JoinPtrLess> JoinSet;

   typedef std::vector<SGSequence*> SequenceSet;

   static const SGPosition NullPos;

public:
   SideGraph();
   ~SideGraph();

   /**
    * Check if join is in the graph and return it else NULL.
    */
   const Join* getJoin(const Join* join) const;

   /**
    * Join added to the graph.  Note that SideGraph will take over
    * responsability for freeing memory of the join.
    * Duplicate joins are silently filtered.  
    */
   const Join* addJoin(const Join* join);

   /** 
    * Get the join set, ie to iterate over or something
    */
   const JoinSet* getJoinSet() const;

   /**
    * Fetch sequence using its ID
    */   
   const SGSequence* getSGSequence(sg_seqid_t id) const;

   /**
    * Get number of sequences
    */
   sg_int_t getNumSequences() const;
   
   /** 
    * Add a sequence into the graph.  SideGraph takes deletion resp 
    * ID field will be ignored and updated automatically
    */
   sg_seqid_t addSGSequence(SGSequence* seq);
   
protected:

   JoinSet _joinSet;
   SequenceSet _seqSet;

private:
   SideGraph(const SideGraph& sg);
};

inline const Join* SideGraph::getJoin(const Join* join) const
{
  JoinSet::const_iterator i = _joinSet.find(join);
  return i != _joinSet.end() ? *i : NULL;
}

inline const Join* SideGraph::addJoin(const Join* join)
{
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

inline const SGSequence* SideGraph::getSGSequence(sg_seqid_t id)
  const
{
  return _seqSet[id];
}

inline sg_int_t SideGraph::getNumSequences() const
{
  return _seqSet.size();
}

inline sg_seqid_t SideGraph::addSGSequence(SGSequence* seq)
{
  seq->_id = _seqSet.size();
  _seqSet.push_back(seq);
  return seq->_id;
}

#endif
