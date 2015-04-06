/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGLOOKUP_H
#define _SGLOOKUP_H

#include <map>

#include "sidegraph.h"

/*
 * Structure to look up a genome coordinate (in its independent
 * coordinate system (ie original Fasta file) as used by Cactus/HAL)
 * in the Side Graph that it's been added to. This structure therefore
 * needs to be updated every time we add to the side graph. 
 */
class SGLookup
{
public:
   SGLookup();
   ~SGLookup();

   /** Set some basic dimensions properties, zapping any existing
    * data.
    */
   void init(const std::vector<std::string>& sequenceNames);

   /** Add a mapping between a pair of intervals:
    * [inPos, inPos+length) -> [outPos, outPos+length) 
    * Note: this is a bijection, and we handle no overlaps.  every
    * source interval maps to exactly one unique target interval. 
    */
   void addInterveral(const SGPosition& inPos,
                      const SGPosition& outPos,
                      sg_int_t length, bool reversed);

   /** Find the position of a genome coordinate in the sequence 
    * graph */
   SGPosition mapPosition(const SGPosition& inPos) const;

   /** Map sequence name (string) to sequence id (integer position in
    * hals sequence list for this genome)
    */
   sg_int_t getSequenceId(const std::string& sequenceName) const;

   /** Map sequence id  (integer position in hals sequence list for this 
    * genome) to its name
    */
   std::string getSequenceName(sg_int_t sequenceId) const;
                 
protected: 

   // we store intervals as sorted list of points.  the leftmost
   // point of the interval is what stores the actual map.  the right
   // point just marks the end of the interval. 
   typedef std::map<sg_int_t, SGPosition> PosMap;
   
   // use array to map sequence id to maybe save some time
   // over storing in single map and having id part of key.
   typedef std::vector<PosMap> PosMapVec;

   // map a sequence name to id
   typedef std::map<std::string, sg_int_t> SeqNameMap;
   // and back
   typedef std::vector<std::string> SeqIdMap;
   
protected:

   PosMapVec _mapVec;
   SeqNameMap _seqNameToId;
   SeqIdMap _seqIdToName;
};

inline sg_int_t SGLookup::getSequenceId(const std::string& sequenceName) const
{
  SeqNameMap::const_iterator i = _seqNameToId.find(sequenceName);
  if (i != _seqNameToId.end())
  {
    return i->second;
  }
  return -1;
}

inline std::string SGLookup::getSequenceName(sg_int_t sequenceId) const
{
  return _seqIdToName.at(sequenceId);
}

inline SGPosition SGLookup::mapPosition(const SGPosition& inPos) const
{
  const PosMap& pm = _mapVec.at(inPos._seqid);
  PosMap::const_iterator i = pm.lower_bound(inPos._pos);
  assert(i != pm.end());

  // todo reversed?
  if (i->first > inPos._pos)
  {
    assert(i != pm.begin());
    --i;
  }
  assert(i->first <= inPos._pos);
  assert(i->second != SideGraph::NullPos);
  
  sg_int_t offset = inPos._pos - i->first;
  SGPosition outPos = i->second;
  outPos._pos += offset;
  return outPos;
}

#endif
