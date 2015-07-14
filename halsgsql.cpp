/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "md5.h"
#include "halsgsql.h"

using namespace std;
using namespace hal;

HALSGSQL::HALSGSQL() : SGSQL(), _sgBuilder(0), _writeAncestralPaths(false)
{
}

HALSGSQL::~HALSGSQL()
{
}

void HALSGSQL::exportGraph(const SGBuilder* sgBuilder,
                           const string& sqlInsertPath,
                           const string& fastaPath, const string& halPath,
                           bool writeAncestralPaths)
{
  _sgBuilder = sgBuilder;
  _halPath = halPath;
  _writeAncestralPaths = writeAncestralPaths;

  writeDb(sgBuilder->getSideGraph(), sqlInsertPath, fastaPath);
}

void HALSGSQL::getSequenceString(const SGSequence* seq,
                                 std::string& outString) const
{
  _sgBuilder->getSequenceString(seq, outString);
}

string HALSGSQL::getOriginName(const SGSequence* seq) const
{
  return _sgBuilder->getHalGenomeName(seq);
}

string HALSGSQL::getPrimaryOriginName() const
{
  return _sgBuilder->getPrimaryGenomeName();
}

string HALSGSQL::getDescription() const
{
  return string("hal2sg ") + _halPath;
}

/*
CREATE TABLE VariantSet (ID INTEGER PRIMARY KEY,
	referenceSetID INTEGER NOT NULL REFERENCES ReferenceSet(ID),
	name TEXT);

CREATE TABLE Allele (ID INTEGER PRIMARY KEY, 
	variantSetID INTEGER REFERENCES VariantSet(ID), 
	name TEXT); -- Naming the allele is optional
--
CREATE TABLE AllelePathItem (alleleID INTEGER REFERENCES allele(ID), 
	pathItemIndex INTEGER NOT NULL, -- zero-based index of this pathItem within the entire path
	sequenceID INTEGER NOT NULL REFERENCES Sequence(ID), 
	start INTEGER NOT NULL,
	length INTEGER NOT NULL, 
	strandIsForward BOOLEAN NOT NULL,
	PRIMARY KEY(alleleID, pathItemIndex));
*/
void HALSGSQL::writePathInserts()
{
  // For every genome in the input HAL, we create one Variant Set.
  // For every sequence in a genome, we create one
  // Allele, and a list of Allele Path Items.

  // TODO: refactor so that formatting logic gets de-coupled from hal
  // and sgbuilder and moved to sgExport/sgsql.cpp...

  vector<const Sequence*> halSequences = _sgBuilder->getHalSequences();
  map<const Genome*, size_t> genomeIdMap;

  // filter out ancestral genomes if not wanted
  if (_writeAncestralPaths == false)
  {
    vector<const Sequence*> leafSequences;
    for (size_t i = 0; i < halSequences.size(); ++i)
    {
      if (halSequences[i]->getGenome()->getNumChildren() == 0)
      {
        leafSequences.push_back(halSequences[i]);
      }
    }
    swap(halSequences, leafSequences);
  }

  // create a variant set for every genome
  for (size_t i = 0; i < halSequences.size(); ++i)
  {
    const Sequence* halSeq = halSequences[i];
    pair<map<const Genome*, size_t>::iterator, bool> ret = genomeIdMap.insert(
      pair<const Genome*, size_t>(halSeq->getGenome(),
                                  genomeIdMap.size()));
    if (ret.second == true)
    {
      _outStream << "INSERT INTO VariantSet VALUES ("
                 << ret.first->second << ", "
                 << 0 << ", "
                 << "'" << halSeq->getGenome()->getName() << "');\n";
    }
  }
  _outStream << endl;

  // create an allele for every sequence; 
  for (size_t i = 0; i < halSequences.size(); ++i)
  {
    _outStream << "INSERT INTO Allele VALUES ("
               << i << ", "
               << genomeIdMap.find(halSequences[i]->getGenome())->second << ", "
               << "'" << _sgBuilder->getHalSeqName(halSequences[i]) << "'" 
               << ");\n";
  }
  _outStream << endl;

  // create a path (AellePathItem) for every sequence
  for (size_t i = 0; i < halSequences.size(); ++i)
  {
    _outStream << "-- PATH for HAL input sequence "
               << _sgBuilder->getHalSeqName(halSequences[i]) << "\n";
    vector<SGSegment> path;

    // note we're calling this a second time here.  if this is costly
    // we can cache the paths or compute the joins once on the fly here
    // (because this was done a first time in computeJoins())
    _sgBuilder->getHalSequencePath(halSequences[i], path);
    for (size_t j = 0; j < path.size(); ++j)
    {
      _outStream << "INSERT INTO AllelePathItem VALUES ("
                 << i << ", "
                 << j << ", "
                 << path[j].getSide().getBase().getSeqID() << ", "
                 << path[j].getSide().getBase().getPos() << ", "
                 << path[j].getLength() << ", "
                 << (path[j].getSide().getForward() ? "\'TRUE\'" : "\'FALSE\'")
                 << ");\n";
      
    }
    _outStream <<endl;
  }

}

