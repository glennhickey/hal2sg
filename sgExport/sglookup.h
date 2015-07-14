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
 * needs to be updated every time we add to the side graph.  Note we only
 * map a single (HAL) genome into a sidegraph here (so need to keep one
 * of these going for each genome added).   
 *
 * Note: Using SGPositions to represent HAL coordinates in interface.  
 * The Sequence ID of a HAL sequence is its array index:
 * Sequence::getArrayIndex().  
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

   /** Allow interface (for sglookback) to add new sequences after init
    */
   void addSequence(const std::string& sequenceName);

   /** Add a mapping between a pair of intervals:
    * [inPos, inPos+length) -> [outPos, outPos+length) 
    * Note: this is a bijection, and we handle no overlaps.  every
    * source interval maps to exactly one unique target interval. 
    */
   void addInterval(const SGPosition& inPos,
                    const SGPosition& outPos,
                    sg_int_t length, bool reversed);

   /** Find the position of a genome coordinate in the sequence 
    * graph.  If outDist is specified, it will return the distance
    * between the output position, and the next waypoint in the 
    * lookup table (how long the interval can extend before potentially
    * spanning a rearrangement on the target) */
   SGSide mapPosition(const SGPosition& inPos,
                      sg_int_t* outDist = NULL,
                      bool outDistReversed = false) const;

   /** Get a path of an inclusive range in a single HAL sequence
    * through the side graph.   
    */
   void getPath(const SGPosition& startPos,
                const SGPosition& endPos,
                std::vector<SGSegment>& outPath) const;
                 
protected: 

   // we store intervals as sorted list of points.  the leftmost
   // point of the interval is what stores the actual map.  the right
   // point just marks the end of the interval. 
   typedef std::map<sg_int_t, SGSide> PosMap;
   
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

   friend std::ostream& operator<<(std::ostream& os, const SGLookup& sg);
};


inline std::ostream& operator<<(std::ostream& os, const SGLookup& sg)
{
  os << "SGLookup:\n";
  for (size_t i = 0; i < sg._mapVec.size(); ++i)
  {
    os << i << ": ";
    for (std::map<sg_int_t, SGSide>::const_iterator j = sg._mapVec[i].begin();
         j != sg._mapVec[i].end(); ++j)
    {
      os << j->first;
      if (j->second.getForward() == false)
      {
        os << "r";
      }
      os << "[" << j->second << "], ";
    }
    os << "\n";
  }
  return os;
}



#endif
