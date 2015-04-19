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
 * The normal lookup structure can map a HAL position onto the side graph. 
 * This lookup structure, further maps the position on a side graph to 
 * "homologous" positions with different bases (ie created by previous)
 * SNP bubbles.  This way we can quickly check, when adding a SNP,
 * if we need to make a new sequence, or if a sequence containing that 
 * base is already there for us to hook into. 
 *
 * We only keep a SNP list at positions where there is a snp.  
 * And for now we use a single SNP handler to index all SNPs in the graph.
 *
 * Also, note that multiple positions in the side graph can correspond
 * to the same SNP alignment.  In this case these positions will all
 * map to the same SNPList pointer. 
 *
 */
class SNPHandler
{
public:

   SNPHandler(bool caseSensitive = false);
   ~SNPHandler();
   
   /** Check to see if SNP present in Side Graph.  If it's not then
    * SideGraph::NullPos is returned */
   const SGPosition& findSNP(const SGPosition& pos, char nuc);

   /** Add a snp to the lookup structure.  This doesn't add it to the graph
    * (ie create necessary sequence and joins -- that has to be done 
    * elsewhere */
   void addSNP(const SGPosition& pos, char nuc,
               const SGPosition& snpPosition);

   /** Get previous SNP position in map in relation to last call to 
    * findSNP.  Saves a lookup if we want to check a series of 
    * contiguous positions, for example.  throws an exception if
    * there is not last position. 
    */
   const SGPosition& getPrevSNP() const;

   /** Get next SNP position in map in relation to last call to 
    * findSNP 
    */
   const SGPosition& getNextSNP() const;

protected:

   struct SNP
   {
      SNP();
      SNP(const SGPosition& pos, char nuc);
      SGPosition _pos;
      char _nuc;
   };

   typedef std::vector<SNP> SNPList;
   typedef std::map<SGPosition, SNPList*> SNPMap;
   
protected:

   bool _caseSens;
   SNPMap _snpMap;
   SNPMap::iterator _cacheIt;
   SGPosition _cachePos;
};


inline SNPHandler::SNP::SNP() {}
inline SNPHandler::SNP::SNP(const SGPosition& pos, char nuc) :
  _pos(pos), _nuc(nuc){}

#endif
