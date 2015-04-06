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
                  hal_index_t start = hal::NULL_INDEX,
                  hal_index_t end = hal::NULL_INDEX);

             
protected:

   /** Find the nearest genome in the Side Graph to align to
    */
   const hal::Genome* getTarget(const hal::Genome* genome);

   /** Map a Sequence onto the Side Graph by aligning to target
    */
   void mapSequence(const hal::Sequence* sequence,
                    hal_index_t globalStart,
                    hal_index_t globalEnd,
                    const hal::Genome* target);

   typedef std::map<std::string, SGLookup*> GenomeLUMap;

protected:
   
   SideGraph* _sg;
   const hal::Genome* _root;
   hal::AlignmentConstPtr _alignment;
   GenomeLUMap _luMap;
   
};
#endif
