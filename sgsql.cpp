/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "sgsql.h"

using namespace std;
using namespace hal;

SGSQL::SGSQL() : _sgBuilder(0), _sg(0)
{
}

SGSQL::~SGSQL()
{
}

void SGSQL::writeDb(const SGBuilder* sgBuilder, const string& sqlInsertPath,
                    const string& fastaPath)
{
  _sgBuilder = sgBuilder;
  _sg = _sgBuilder->getSideGraph();
  assert(_sg != NULL);
  _outPath = sqlInsertPath;
  _outStream.open(_outPath);
  assert(_outStream);
  _faPath = fastaPath;
  _faStream.open(_faPath);
  assert(_faStream);
  _checksumMap.clear();

  writeFasta();
  writeFastaInsert();
  writeSequenceInserts();
  writeJoinInserts();
  writePathInserts();

  _outStream.close();
  _faStream.close();
}

void SGSQL::writeFasta()
{
  string dnaBuffer;
  for (sg_int_t i = 0; i < _sg->getNumSequences(); ++i)
  {
    const SGSequence* seq = _sg->getSequence(i);
    _sgBuilder->getSequenceString(seq, dnaBuffer);
    assert(dnaBuffer.length() == seq->getLength());
    _faStream << ">"  << seq->getName() << "\n";
    size_t len;
    for (size_t pos = 0; pos < dnaBuffer.length(); pos += len)
    {
      len = min((size_t)80, dnaBuffer.length() - pos);
      _faStream << dnaBuffer.substr(pos, len) << "\n";
    }
    string checksum;
    getChecksum(dnaBuffer, checksum);
    _checksumMap.insert(pair<sg_seqid_t, string>(seq->getID(), checksum));
  }
}

/*
  CREATE TABLE FASTA (ID INTEGER PRIMARY KEY,
  fastaURI TEXT NOT NULL);
*/
void SGSQL::writeFastaInsert()
{
  _outStream << "INSERT INTO FASTA VALUES ("
             << 0 << ", "
             << "file://" << _faPath << ")\n";
  _outStream << endl;
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
    map<sg_seqid_t, string>::const_iterator j = _checksumMap.find(seq->getID());
    assert(j != _checksumMap.end());
    _outStream << "INSERT INTO Sequence VALUES ("
               << seq->getID() << ", "
               << 0 << ", "
               << "\'"<< seq->getName() << "\', "
               << "\'" << j->second  << "\', "
               << seq->getLength() << ");\n";
  }
  _outStream << endl;
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
    _outStream << "INSERT INTO GraphJoin VALUES ("
               << count++ << ", "
               << side1.getBase().getSeqID() << ", "
               << side1.getBase().getPos() << ", "
               << (side1.getForward() ? "TRUE" : "FALSE") << ", "
               << side2.getBase().getSeqID() << ", "
               << side2.getBase().getPos() << ", "
               << (side2.getForward() ? "TRUE" : "FALSE") << ")\n";
  }
  _outStream << endl;
}


void SGSQL::writePathInserts()
{
  const SGBuilder::PathMap* pathMap = _sgBuilder->getPathMap();
  SGBuilder::PathMap::const_iterator i;
  _outStream << "# (HAL_Genome.SequenceName, Step#, SequenceID, Position, "
             << "Forward)" << endl;
  for (i = pathMap->begin(); i != pathMap->end(); ++i)
  {
    const Sequence* halSeq = i->first;
    const SGBuilder::SidePath* path = i->second;
    _sgBuilder->verifyPath(halSeq, path);
    SGBuilder::SidePath::const_iterator j;
    for (size_t j = 0; j < path->size(); ++j)
    {
      const SGSide& side = path->at(j);
      _outStream << "INSERT INTO Path VALUES ("
                 << halSeq->getFullName() << ", "
                 << j << ", "
                 << side.getBase().getSeqID() << ", "
                 << side.getBase().getPos() << ", "
                 << (side.getForward() ? "TRUE" : "FALSE") << ")\n";
    }
  }
  _outStream <<endl;
}

void SGSQL::getChecksum(const string& inputString, string& outputString)
{
  outputString = "md5";
}
