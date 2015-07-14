/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _HALSGSQL_H
#define _HALSGSQL_H

#include <cstdlib>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

#include "sgsql.h"
#include "sgbuilder.h"


/*
 * write a SideGraph to GA4GH SQL format.  This will be a fasta file with
 * the sequence information and a .sql file with a bunch of INSERT statements
 * that can be loaded into a db.  This class has the HAL conversion 
 * specific logic coming from SGBuilder...
 */

class HALSGSQL : public SGSQL
{
public:
   HALSGSQL();
   virtual ~HALSGSQL();

   /** write out the graph as a database 
    */
   void exportGraph(const SGBuilder* sgBuilder,
                    const std::string& sqlInsertPath,
                    const std::string& fastaPath, const std::string& halPath,
                    bool writeAncestralPaths = true);
protected:


   /** write path INSERTs (makes a VariantSet for each Genome and 
    * an Allele for each sequence
    */
   void writePathInserts();

   /** get DNA string corresponding to a sequence 
    */
   void getSequenceString(const SGSequence* seq,
                          std::string& outString) const;

   /** determine name of input sequence from which this sequence was
    * derived.  
    */
   std::string getOriginName(const SGSequence* seq) const;

   /** Which origin name do we consider primary? Used to name ReferenceSet
    */
   std::string getPrimaryOriginName() const;

   /** Get description for reference set 
    */
   std::string getDescription() const;
   
protected:

   const SGBuilder* _sgBuilder;
   bool _writeAncestralPaths;
};


#endif
