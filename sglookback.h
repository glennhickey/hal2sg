/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGLOOKBACK_H
#define _SGLOOKBACK_H

#include "hal.h"
#include "sglookup.h"

/*
 * Wrap SGLookup structure to make something that could map in other
 * direction.  Issue is that hal sequence ids are not globally unique 
 * (rather only unique per genome).  So we add a little table to keep track
 * of this rather than changing SGLookup itself.. The other change is not
 * being able to specify number of sequences in init. 
 */
class SGLookBack
{
public:
   SGLookBack();
   ~SGLookBack();

   void clear();
   /** Add a mapping between a SideGraph interval and a HAL interval
    * Note: this is a bijection, and we handle no overlaps.  every
    * source interval maps to exactly one unique target interval. 
    */
   void addInterval(const SGPosition& inPos, const hal::Sequence* outSeq,
                    hal_index_t outOffset, 
                    sg_int_t length, bool reversed);

   /** Find the position of a SideGraph coordinate in the HAL graph.  */
   std::pair<const hal::Sequence*, SGSide> mapPosition(
     const SGPosition& inPos) const;

   /** Get a path of an inclusive range in a single SG sequence
    * through the HAL graph.   
    */
   void getPath(const SGPosition& startPos,
                const SGPosition& endPos,
                std::vector<SGSegment>& outPath,
                std::vector<const hal::Sequence*>& outHalSeqs) const;

   std::string getHalGenomeName(const SGSequence* sgSeq) const;
   
protected:

   std::vector<const hal::Sequence*> _halSeqMap;
   SGLookup _lookup;
};

inline std::string SGLookBack::getHalGenomeName(const SGSequence* sgSeq) const
{
  assert(sgSeq->getID() < _halSeqMap.size());
  return _halSeqMap[sgSeq->getID()]->getGenome()->getName();
}

#endif
