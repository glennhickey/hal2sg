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
   void init(hal::AlignmentConstPtr alignment, const hal::Genome* root = NULL,
             bool referenceDupes = true, bool noSubMode = false);

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

   /**
    * Get DNA bases for a Side Graph sequence by querying the
    * source sequence in the HAL file. 
    */
   size_t getSequenceString(const SGSequence* sgSequence,
                            std::string& outString,
                            sg_int_t pos = 0,
                            sg_int_t length = -1) const;
   

   /** Returns true if sequence was created from the first genome added
    * (root in current impl), true otherwise */
   const std::string& getPrimaryGenomeName() const;

   /** Look up the HAL genome that he Side Graph sequence was based 
    */
   std::string getHalGenomeName(const SGSequence* sgSequence) const;
   

   typedef std::vector<SGSide> SidePath;
   struct SeqLess {
      bool operator()(const hal::Sequence* s1, const hal::Sequence* s2) const;
   };
   typedef std::map<const hal::Sequence*, SidePath*, SGBuilder::SeqLess>
   PathMap;

   /** Get the paths
    */
   const PathMap* getPathMap() const;

   /** Check to make sure the path exactly covers the sequence and 
    * all the joins exist.  Assertion triggered if something's amiss
    * (and returns false).  May not be the fastest function.  */
   bool verifyPath(const hal::Sequence* sequence, const SidePath* path) const;
       

protected:
   
   typedef std::map<const SGSequence*, std::pair<const hal::Sequence*,
                                                 hal_index_t> > SequenceMapBack;

   typedef std::map<std::string, SGLookup*> GenomeLUMap;

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
   const SGJoin* createSGJoin(const SGSide& side1, const SGSide& side2);

   /** Add a step to the path */
   void addPathStep(const SGSide& side);

   /** Add interval (from blockmapper machinery) to the side graph.  
    * The interval maps from the new SOURCE genome to a TARGET genome
    * that is already in the side graph. 
    */
   void processBlock(Block* prevBlock,
                     Block* block,
                     Block* nextBlock,
                     SGSide& prevHook,
                     const hal::Sequence* srcSequence,
                     const hal::Genome* srcGenome,
                     hal_index_t sequenceStart,
                     hal_index_t sequenceEnd,
                     const hal::Genome* tgtGenome);

   /** Add a block, breaking apart for SNPs. only add joins that are
    * contained in the block.  Update the endpoints if needed (as result of 
    * snps) */
   void mapBlockSnps(const Block*, std::pair<SGSide, SGSide>& blockEnds);

   /** Add a slice of a block.  Either every base is a snp (snp==true) or
    * no bases are a snp.  hook on the previous hook. */
   std::pair<SGSide, SGSide> mapSliceSnps(const Block* block,
                                          hal_index_t srcStartOffset,
                                          hal_index_t srcEndOffset,
                                          bool snp, SGSide& hook,
                                          bool sgForwardMap,
                                          const std::string& srcDNA,
                                          const std::string& tgtDNA);

   /**
    * Don't want to deal with the mapped block fragments all the time. 
    * code to read into struct here. 
    */
   void fragmentsToBlock(
     const std::vector<hal::MappedSegmentConstPtr>& fragments,
     Block& block) const;

   /** check if block aligns something to itself */
   bool isSelfBlock(const Block& block) const;

   /** cut block against another */
   Block* cutBlock(Block* prev, Block* cur);

   /** We are anchoring on the root genome (at least for now).  But in
    * Adams output, the root sequence is Ns which is a problem.  We 
    * use the function as an override to map a root sequence from its children
    * using the knowledge that there are no substitutions */
   void getRootSubString(std::string& outDNA, const hal::Sequence* sequence,
                         hal_index_t pos, hal_index_t length) const;

protected:

   SideGraph* _sg;
   const hal::Genome* _root;
   hal::AlignmentConstPtr _alignment;
   GenomeLUMap _luMap;
   SGLookup* _lookup;
   SequenceMapBack _seqMapBack;
   std::set<const hal::Genome*> _mapPath;
   const hal::Genome* _mapMrca;
   bool _referenceDupes;
   SGJoin* _lastJoin;
   SidePath* _path;
   PathMap _pathMap;
   bool _inferRootSeq;
   mutable std::string _rootString;
   bool _noSubMode;
   size_t _pathLength;
   size_t _joinPathLength;
   size_t _sgJoinPathLength;
   std::string _firstGenomeName;

   friend std::ostream& operator<<(std::ostream& os, const Block* block);

};

std::ostream& operator<<(std::ostream& os, const SGBuilder::Block* block);

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

inline std::ostream& operator<<(std::ostream& os, const SGBuilder::Block* block)
{
  if (block == NULL)
  {
    os << "BLOCK(NULL)";
  }
  else
  {
    os << "BLOCK(" << block->_srcSeq->getName() << ": "
       << block->_srcStart << "," << block->_srcEnd <<  ") TGT("
       << block->_tgtSeq->getName() << ": "
       << block->_tgtStart << "," << block->_tgtEnd << ") "
       << block->_reversed;
  }
  return os;
}

inline bool SGBuilder::isSelfBlock(const SGBuilder::Block& block) const
{
  return block._srcStart == block._tgtStart &&
     block._srcEnd == block._tgtEnd &&
     block._srcSeq == block._tgtSeq;
}

inline bool SGBuilder::SeqLess::operator()(const hal::Sequence* s1,
                                           const hal::Sequence* s2) const
{
  return s1->getFullName() < s2->getFullName();
}

inline const SGBuilder::PathMap* SGBuilder::getPathMap() const
{
  return &_pathMap;
}


#endif
