/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "md5.h"
#include "sgsql.h"

using namespace std;


SGSQL::SGSQL() :  _sg(0)
{
}

SGSQL::~SGSQL()
{
}

void SGSQL::writeDb(const SideGraph* sg, const string& sqlInsertPath,
                    const string& fastaPath)
{
  _sg = sg;
  assert(_sg != NULL);
  _outPath = sqlInsertPath;
  _outStream.open(_outPath.c_str());
  assert(_outStream);
  _faPath = fastaPath;
  _faStream.open(_faPath.c_str());
  assert(_faStream);

  _checksumMap.clear();

  writeFasta();
  writeFastaInsert();
  writeSequenceInserts();
  writeReferenceInserts();
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
    getSequenceString(seq, dnaBuffer);
    assert((sg_int_t)dnaBuffer.length() == seq->getLength());
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
             << "'" << _faPath << "');\n";
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
-- References
CREATE TABLE Reference (ID INTEGER PRIMARY KEY, 
	name TEXT NOT NULL,
	updateTime DATE NOT NULL, 
	sequenceID INTEGER NOT NULL,
	start INTEGER NOT NULL, -- default 0.
	length INTEGER, -- if null, assume sequence.lenght - start 
	md5checksum TEXT, -- if null, assume sequence.md5checksum
	isDerived BOOLEAN NOT NULL, 
	sourceDivergence REAL, -- may be null if isDerived is FALSE (not sure if it needs to be stored)?
	ncbiTaxonID INTEGER, -- ID from http://www.ncbi.nlm.nih.gov/taxonomy, may be null
	isPrimary BOOLEAN NOT NULL,
	FOREIGN KEY(sequenceID) REFERENCES Sequence(ID)); 
CREATE TABLE ReferenceAccession (ID INTEGER PRIMARY KEY,
	referenceID INTEGER NOT NULL,
	accessionID TEXT NOT NULL,
	FOREIGN KEY(referenceID) REFERENCES Reference(ID)); 

-- Reference sets
CREATE TABLE ReferenceSet (ID INTEGER PRIMARY KEY,
	ncbiTaxonID INT, -- may be null, may differ from ncbiTaxonID of contained Reference record
	description TEXT, -- may be null, but shouldn't be?
	fastaID INTEGER, -- may be null. TODO: What if both this and a member Reference's Sequence fastaID are null?
	assemblyID TEXT, -- may be null, but REALLY shouldn't be?
	isDerived BOOLEAN NOT NULL,
	FOREIGN KEY(fastaID) REFERENCES FASTA(ID));
CREATE TABLE ReferenceSetAccession (ID INTEGER PRIMARY KEY,
	referenceSetID INTEGER NOT NULL,
	accessionID TEXT NOT NULL,
	FOREIGN KEY(referenceSetID) REFERENCES ReferenceSet(ID)); 

CREATE TABLE Reference_ReferenceSet_Join (referenceID INTEGER NOT NULL, 
	referenceSetID INTEGER NOT NULL,
	PRIMARY KEY(referenceID, referenceSetID),
	FOREIGN KEY(referenceID) REFERENCES Reference(ID),
	FOREIGN KEY(referenceSetID) REFERENCES ReferenceSet(ID));
*/
void SGSQL::writeReferenceInserts()
{
  // make a single reference set
  _outStream << "INSERT INTO ReferenceSet VALUES "
             << "(0, NULL, "
             << "\'" << getDescription() << "\'"
             << ", \'" << getPrimaryOriginName() << "\'"
             << ", \'FALSE\');\n";
  
  for (sg_int_t i = 0; i < _sg->getNumSequences(); ++i)
  {
    const SGSequence* seq = _sg->getSequence(i);
    string halGenomeName = getOriginName(seq);
    bool primary = halGenomeName == getPrimaryOriginName();
    _outStream << "INSERT INTO Reference VALUES ("
               << seq->getID() << ", "
               << "\'"<< seq->getName() << "\', "
               << "date(\'now\')" << ", "
               << seq->getID() << ", "
               << 0 << ", "
               << "NULL " << ", "
               << "NULL " << ", "
               << "\'FALSE\'" << ", "
               << "NULL " << ", "
               << "NULL " << ", "
               << (primary ? "\'TRUE\'" : "\'FALSE\'")
               << ");\n";
  }

  for (sg_int_t i = 0; i < _sg->getNumSequences(); ++i)
  {
    const SGSequence* seq = _sg->getSequence(i);
    _outStream << "INSERT INTO Reference_ReferenceSet_Join VALUES ("
               << seq->getID() << ", 0);\n";
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
               << (side1.getForward() ? "'TRUE'" : "'FALSE'") << ", "
               << side2.getBase().getSeqID() << ", "
               << side2.getBase().getPos() << ", "
               << (side2.getForward() ? "'TRUE'" : "'FALSE'") << ");\n";
  }
  _outStream << endl;
}



void SGSQL::getChecksum(const string& inputString, string& outputString)
{
  outputString = md5(inputString);
}

