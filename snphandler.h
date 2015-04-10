/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SNPHANDLER_H
#define _SNPHANDLER_H

#include <map>

#include "sglookup.h"

/**
 * Structure to link a position in a sidegraph with alternate bases
 * ie to represent point mutations in the hal.  These mutations
 * need to become snp bubbles in the Side Graph.  But we keep track of
 * them separately as we build.
 *
 * So say we have a sequence AAACCGTT in the sidegraph. 
 * The SNP handler maps each of these positions to a list of SNPs:
 *
 *    C
 *    |
 *   GG
 *   ||
 *   CT  A
 *   ||  |
 *  AAACCGTT
 *
 * We always need to be able to walk columns in this structure to make sure
 * that the same value doesn't appear more than once per column. 
 * When we iterate over a HAL genome, we effectively drop new bases onto
 * the top of this structure.  We have enough information to iteratively 
 * add sequences and joins as necessary to keep values unique and make sure
 * that multibase snps get merged into one. 
 *
 */
class SNPHandler
{
public:

   SNPHandler();
   ~SNPHandler();

   void addSNP(const std::vector<char>
   
public:
   
   struct SNP
   {
      char _nuc;
      SGPosition _pos;
   };
   
protected:

   typedef std::vector<SNP> SNPList;
   typedef std::map<SGPosition, SNPList*> SNPMap;
   
protected:

   SNPMap _snpMap;
   
};

#endif
