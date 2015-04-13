/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGSQL_H
#define _SGSQL_H

#include <cstdlib>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

#include "sidegraph.h"
#include "sgbuilder.h"

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

   /** write out the graph as a database 
    */
   void writeDb(const SGBuilder* sgBuilder, const std::string& sqlInsertPath,
                const std::string& fastaStream);
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

   /** write path INSERTs in my own invented format
    */
   void writePathInserts();

   static void getChecksum(const std::string& inputString,
                           std::string& outputString);

protected:

   const SGBuilder* _sgBuilder;
   const SideGraph* _sg;
   std::string _outPath;
   std::string _faPath;
   std::ofstream _outStream;
   std::ofstream _faStream;
   std::map<sg_seqid_t, std::string> _checksumMap;
};


#endif
