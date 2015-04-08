/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGBUILDER_H
#define _SGBUILDER_H

#include "hal.h"
#include "sidegraph.h"
#include "sglookup.h"

/*
 * Iteratively add genomes from a HAL file to our SideGraph. 
 */
class SGBuilder
{
public:
   SGBuilder();
   ~SGBuilder();

   /** 
    * Set the alignment
    */
   void init(hal::AlignmentConstPtr alignment, const hal::Genome* root = NULL);

   /**
    * Erase everything
    */
   void clear();

   /**
    * Add Genome to the Side Graph.  Can optionally add only a sub-region
    */
   void addGenome(const hal::Genome* genome,
                  const hal::Sequence* sequence = NULL,
                  hal_index_t start = 0,
                  hal_index_t length = 0);

   /**
    * Get the Side Graph
    */
   const SideGraph* getSideGraph() const;

protected:
   
   typedef std::map<std::string, SGLookup*> GenomeLUMap;

   typedef std::map<SGSequence*, std::pair<const hal::Sequence*,
                                           hal_index_t> > SequenceMapBack;

   /** convenience structure for alignment block.  note hall coordinates
    * are in forward strand relative to Segment (not genome).  */
   struct Block {
      hal_index_t _srcStart;
      hal_index_t _srcEnd;
      hal_index_t _tgtStart;
      hal_index_t _tgtEnd;
      const hal::Sequence* _srcSeq;
      const hal::Sequence* _tgtSeq;
      bool _reversed;
   };
   struct BlockPtrLess {
      bool operator()(const Block* b1, const Block* b2) const;
   };

protected:

   /** Find the nearest genome in the Side Graph to align to */
   const hal::Genome* getTarget(const hal::Genome* genome);

   /** Map a Sequence onto the Side Graph by aligning to target
    */
   void mapSequence(const hal::Sequence* sequence,
                    hal_index_t globalStart,
                    hal_index_t globalEnd,
                    const hal::Genome* target);

   /** Add a sequence (or part thereof to the sidegraph) and update
    * lookup structures (but not joins) */
   SGSequence* createSGSequence(const hal::Sequence* sequence,
                                hal_index_t startOffset,
                                hal_index_t length);

   /** Add a join to the side graph */
   const SGJoin* createSGJoin(sg_seqid_t seqId1, sg_int_t pos1, bool forward1,
                              sg_seqid_t seqId2, sg_int_t pos2, bool forward2);

   /** Add interval (from blockmapper machinery) to the side graph.  
    * The interval maps from the new SOURCE genome to a TARGET genome
    * that is already in the side graph. 
    */
   void updateSegment(Block* prevBlock,
                      Block* block,
                      Block* nextBlock,
                      SGPosition* prevHook,
                      const hal::Sequence* srcSequence,
                      const hal::Genome* srcGenome,
                      hal_index_t globalStart, hal_index_t globalEnd,
                      const hal::Genome* tgtGenome);

   /**
    * Don't want to deal with the mapped block fragments all the time. 
    * code to read into struct here. 
    */
   void fragmentsToBlock(
     const std::vector<hal::MappedSegmentConstPtr>& fragments,
     Block& block) const;

protected:
   
   SideGraph* _sg;
   const hal::Genome* _root;
   hal::AlignmentConstPtr _alignment;
   GenomeLUMap _luMap;
   SGLookup* _lookup;
   SequenceMapBack _seqMapBack;
   std::set<const hal::Genome*> _mapPath;
   const hal::Genome* _mapMrca;
   
};

inline bool SGBuilder::BlockPtrLess::operator()(const SGBuilder::Block* b1,
                                                const SGBuilder::Block* b2)
  const
{
  assert(b1->_srcSeq == b2->_srcSeq);
  if (b1->_srcStart < b2->_srcStart)
  {
    return true;
  }
  else if (b1->_srcStart == b2->_srcStart)
  {
    if (b1->_srcEnd < b2->_srcEnd)
    {
      return true;
    }
    else if (b1->_srcEnd == b2->_srcEnd)
    {
      if (b1->_tgtStart < b2->_tgtStart)
      {
        return true;
      }
      else if (b1->_tgtStart == b2->_tgtStart)
      {
        if (b1->_tgtEnd < b2->_tgtEnd)
        {
          return true;
        }
      }
    }
  }
  return false;
}

#endif
