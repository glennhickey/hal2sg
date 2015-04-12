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
#include <cassert>

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

   /** write out the SQL INSERT statements
    */
   void writeInserts(const SideGraph* sideGraph, std::ostream& os);
   
   void writeFasta();
   
protected:

   /** write an INSERT for each join the graph
    */
   void writeJoinInserts();

   /** write an INSERT for each sequence in the graph
    */
   void writeSequenceInserts();

   void writePathInserts();

protected:

   std::ostream* _outStream;
   const SideGraph* _sg;
};


#endif
