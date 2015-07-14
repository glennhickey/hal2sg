/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGSQL_H
#define _SGSQL_H

#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>

#include "sidegraph.h"

//#include <mysql.h>
//#include <mysqld_error.h>

/*
 * write a SideGraph to an SQL server. or, for now, as some insert 
 * commands to a text file. 
 */

class SGSQL
{
public:
   SGSQL();
   virtual ~SGSQL();

protected:

   /** write out the FASTA file by using the back map in sgBuilder
    * to convert sideGraph sequences back into their hal coordinates, then
    * pulling the DNA string out of HAL.  also fill in the checksum map
    * as we go.
    */
   void writeFasta();

   /** write an INSERT for the fasta file.  We only make one so its 
    * ID is always 0
    */
   void writeFastaInsert();
   
   /** write an INSERT for each join the graph
    */
   void writeJoinInserts();

   /** write an INSERT for each sequence in the graph
    */
   void writeSequenceInserts();

   /** write a "Reference" INSERT for each sequence in the graph into
    * one Reference set. 
    */
   void writeReferenceInserts();

   /** write path INSERTs (makes a VariantSet for each Genome and 
    * an Allele for each sequence
    */
   void writePathInserts();

   static void getChecksum(const std::string& inputString,
                           std::string& outputString);

protected:

   const SideGraph* _sg;
   std::string _outPath;
   std::string _faPath;
   std::string _halPath;
   std::ofstream _outStream;
   std::ofstream _faStream;
   std::map<sg_seqid_t, std::string> _checksumMap;
};


#endif
