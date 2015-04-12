/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sgsql.h"

using namespace std;

SGSQL::SGSQL() : _outStream(NULL), _sg(NULL)
{
}

SGSQL::~SGSQL()
{

}

void SGSQL::writeInserts(const SideGraph* sideGraph, ostream& os)
{
  _outStream = &os;
  _sg = sideGraph;

  writeSequenceInserts();
  writeJoinInserts();
}

void SGSQL::writeFasta()
{
}

/*
CREATE TABLE Sequence (ID INTEGER PRIMARY KEY,
	-- FASTA file URI and record name within it containing sequence  
	fastaID INTEGER, -- If null, query enclosing ReferenceSet's fastaID
	sequenceRecordName TEXT NOT NULL, 
	md5checksum TEXT NOT NULL,
	length INTEGER NOT NULL,
	FOREIGN KEY(fastaID) REFERENCES FASTA(ID));
*/
void SGSQL::writeSequenceInserts()
{
  for (sg_int_t i = 0; i < _sg->getNumSequences(); ++i)
  {
    const SGSequence* seq = _sg->getSequence(i);
    *_outStream << "INSERT INTO Sequence VALUES ("
                << seq->getID() << ", "
                << "\'NULL\'" << ", "
                << "\'"<< seq->getName() << "\', "
                << -1 << ", "
                << seq->getLength() << ");\n";
  }
}

/*
CREATE TABLE GraphJoin (ID INTEGER PRIMARY KEY,
	-- startJoin defines parent sequence & position that 5' end attaches to.
	side1SequenceID INTEGER NOT NULL,
	side1Position INTEGER NOT NULL,
	side1StrandIsForward BOOLEAN NOT NULL,
	-- endJoin defines parent sequence & position that 3' end attaches to.
	side2SequenceID INTEGER NOT NULL,
	side2Position INTEGER NOT NULL,
	side2StrandIsForward BOOLEAN NOT NULL,
	FOREIGN KEY(side1SequenceID) REFERENCES Sequence(ID),
	FOREIGN KEY(side2SequenceID) REFERENCES Sequence(ID)); 
*/
void SGSQL::writeJoinInserts()
{
  const SideGraph::JoinSet* joinSet = _sg->getJoinSet();
  SideGraph::JoinSet::const_iterator i;
  size_t count = 0;
  for (i = joinSet->begin(); i != joinSet->end(); ++i)
  {
    const SGJoin* join = *i;
    const SGSide& side1 = join->getSide1();
    const SGSide& side2 = join->getSide2();
    *_outStream << "INSERT INTO GraphJoin VALUES ("
                << count++ << ", "
                << side1.getBase().getSeqID() << ", "
                << side1.getBase().getPos() << ", "
                << (side1.getForward() ? "TRUE" : "FALSE") << ", "
                << side2.getBase().getSeqID() << ", "
                << side2.getBase().getPos() << ", "
                << (side2.getForward() ? "TRUE" : "FALSE") << ")\n";
  }
}

