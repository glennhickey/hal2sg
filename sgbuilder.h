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
#include "sglookback.h"

class SNPHandler;


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
             bool referenceDupes = true, bool camelMode = false,
             bool onlySequenceNames = false);

   /**
    * Erase everything
    */
   void clear();

   /**
    * Erase everything, except the sidegraph pointer is returned but not 
    * deleted (which is responsability of client)
    */
   SideGraph* clear_except_sg();

   /**
    * Add Genome to the Side Graph. sequence, start and length parameters
    * should always be their defaults for now. 
    */
   void addGenome(const hal::Genome* genome,
                  const hal::Sequence* sequence = NULL,
                  hal_index_t start = 0,
                  hal_index_t length = 0);

   /**
    * Joins are computed in a second pass, after all genomes have been
    * added.  This pass will also performa a sanity check to make sure
    * the graph contains a path for each input sequence */
   void computeJoins(bool doAncestralJoins = true);

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
   
   /** Get a list of the hal sequences in order they were processed 
    */
   const std::vector<const hal::Sequence*>& getHalSequences() const;

   /** Get a Segment Path for a given halSequence 
    */
   void getHalSequencePath(const hal::Sequence* halSeq,
                           std::vector<SGSegment>& outPath) const;

   const std::string getHalSeqName(const hal::Sequence* halSeq) const;
   
public:

   struct SeqLess {
      bool operator()(const hal::Sequence* s1, const hal::Sequence* s2) const;
   };
   
protected:
   
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

   /** Compute the alignment blocks between a (sub)sequence and a 
    * target genome */
   void computeBlocks(const hal::Sequence* sequence,
                      hal_index_t globalStart,
                      hal_index_t globalEnd,
                      const hal::Genome* target,
                      std::vector<Block*>& blocks);

   /** Add a sequence (or part thereof to the sidegraph) and update
    * lookup structures (but not joins) */
   std::pair<SGSide, SGSide> createSGSequence(const hal::Sequence* sequence,
                                              hal_index_t startOffset,
                                              hal_index_t length);
   
   /** Add a join to the side graph */
   const SGJoin* createSGJoin(const SGSide& side1, const SGSide& side2);

   /** New function added to wrap mapBlockEnds, which can now be used
    * for self-alignments and normal alignments with slightly different
    * logic */
   void visitBlock(Block* prevBlock,
                   Block* block,
                   Block* nextBlock,
                   SGSide& prevHook,
                   const hal::Sequence* srcSequence,
                   const hal::Genome* srcGenome,
                   hal_index_t sequenceStart,
                   hal_index_t sequenceEnd,
                   const hal::Genome* tgtGenome);
   
   /** Add interval (from blockmapper machinery) to the side graph.  
    * The interval maps from the new SOURCE genome to a TARGET genome
    * that is already in the side graph. 
    */
   std::pair<SGSide, SGSide> mapBlockEnds(const Block* block);

   /** Add a block, breaking apart for SNPs. only add joins that are
    * contained in the block.  Update the endpoints if needed (as result of 
    * snps) */
   std::pair<SGSide, SGSide>
   mapBlockBody(const Block*, const std::pair<SGSide, SGSide>& sgBlockEnds);

   /** Add a slice of a block.  Either every base is a snp (snp==true) or
    * no bases are a snp.  hook on the previous hook. */
   std::pair<SGSide, SGSide>
   mapBlockSlice(const Block* block,
                 const std::pair<SGSide, SGSide>& sgBlockEnds,
                 hal_index_t srcStartOffset,
                 hal_index_t srcEndOffset,
                 bool snp,
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

   /** cut sorted blocks list so blocks either dont overlap or completely
    * overlap (all based on src coordinates) */
   void cutBlocks(std::vector<Block*>&, bool leaveExactOverlaps = false);

   /** Add joins (and do sanity check) for one path corresponding to
    * one input hal sequence */
   void addPathJoins(const hal::Sequence* sequence,
                     const std::vector<SGSegment>& path);

   /** We are anchoring on the root genome (at least for now).  But in
    * Adams output, the root sequence is Ns which is a problem.  We 
    * use the function as an override to map a root sequence from its 
    * children using the knowledge that there are no substitutions */
   void getRootSubString(std::string& outDNA, const hal::Sequence* sequence,
                         hal_index_t pos, hal_index_t length) const;

   
   // Logic added only for collapsing new sequences (ie handling dupes
   // in reference or insertions):
   

   /** When computing duplications on a self-alignment, we want to 
    * pick out blocks that do not get collapsed out due to alignment,
    * as they will be present in the new sequence */
   void getCollapsedFlags(const std::vector<Block*>& blocks,
                          std::vector<bool>& collapseBlock);

   /** Add gaps and uncollapsed blocks to the lookup structure to 
    * make sure the entire new sequence is covered */
   void updateDupeBlockLookups(const hal::Sequence* sequence,
                               hal_index_t startOffset,
                               hal_index_t length,
                               const std::vector<Block*>& blocks,
                               const std::vector<bool>& collapsed,
                               const SGSequence* sgSeq);

   /** In-place filtering of self-alignment blocks that are not needed
    * to construct the graph as their alignment is already covered 
    * by some other block.  Logic used is that blocks with targets
    * not in the lookup are filtered */
   void filterRedundantDupeBlocks(std::vector<Block*>& blocks,
                                  const SGSequence* sgSeq);
   
protected:

   SideGraph* _sg;
   const hal::Genome* _root;
   const hal::Genome* _mapRoot;
   hal::AlignmentConstPtr _alignment;
   GenomeLUMap _luMap;
   SGLookup* _lookup;
   SGLookBack _lookBack;
   std::set<const hal::Genome*> _mapPath;
   const hal::Genome* _mapMrca;
   bool _referenceDupes;
   bool _inferRootSeq;
   mutable std::string _rootString;
   bool _camelMode;
   size_t _pathLength;
   std::string _firstGenomeName;
   std::vector<const hal::Sequence*> _halSequences;
   SNPHandler* _snpHandler;
   bool _onlySequenceNames;

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

inline std::ostream& operator<<(std::ostream& os,
                                const SGBuilder::Block* block)
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
       << "len=" << (block->_srcEnd - block->_srcStart + 1)
       << " rev=" << block->_reversed;
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

inline const std::string SGBuilder::getHalSeqName(const hal::Sequence*
                                                  halSeq) const
{
  return _onlySequenceNames ? halSeq->getName() : halSeq->getFullName();
}



#endif
