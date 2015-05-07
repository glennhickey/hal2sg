/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <limits>
#include <cctype>
#include <cstdio>
#include <sstream>

#include "snphandler.h"

using namespace std;

SNPHandler::SNPHandler(SideGraph* sideGraph, bool caseSensitive)
  :  _caseSens(caseSensitive), _sg(sideGraph)
{

}

SNPHandler::~SNPHandler()
{
  // simple way of making sure we don't delete same pointer twice:
  // add to a set first.
  set<SNPList*> snpSet;
  for (SNPMap::iterator i = _snpMap.begin(); i != _snpMap.end(); ++i)
  {
    snpSet.insert(i->second);
  }
  for (set<SNPList*>::iterator i = snpSet.begin(); i != snpSet.end(); ++i)
  {
    delete *i;
  }
}

pair<SGSide, SGSide> SNPHandler::createSNP(const string& dnaString,
                                           size_t dnaOffset,
                                           size_t dnaLength,
                                           const SGPosition& srcPos,
                                           const SGPosition& tgtPos,
                                           bool reverseMap,
                                           SGLookup* srcLookup)
{  
  vector<SGPosition> sgPositions(dnaLength);
  // note : should use hints / cache to speed up consecutive bases
  // but start simple for baseline tests

  // first we use the snp structure to find out if any of the SNPs we
  // want to add already exist.  If they do, we update sgPositions
  SGPosition tgt(tgtPos);
  for (sg_int_t i = 0; i < dnaLength; ++i)
  {
    if (reverseMap == false)
    {
      tgt.setPos(tgtPos.getPos() + dnaOffset + i);
    }
    else
    {
      tgt.setPos(tgtPos.getPos() - dnaOffset - i);
    }
    sgPositions[i] = findSNP(tgt, dnaString[i + dnaOffset]);
  }

  // now we have a position for all existing snps in the graph.  any
  // positions that are "null" means we need to add new sequences to
  // cover them.  We create the sequences and joins here and add them
  // to the side graph.  We also update sgPositions with these new
  // coordinates.
  SGSide prevHook;
  string nameBuf;
  
  for (sg_int_t i = 0; i < dnaLength; ++i)
  {
    if (sgPositions[i] == SideGraph::NullPos)
    {
      sg_int_t j = i + 1;
      for (; j < dnaLength && sgPositions[j] == SideGraph::NullPos; ++j);
      --j;
      getSNPName(tgtPos, dnaString, dnaOffset, reverseMap, i, nameBuf);
      SGSequence* newSeq = new SGSequence(-1, j - i + 1, nameBuf);
      _sg->addSequence(newSeq);

      if (i > 0)
      {
        SGJoin* inJoin = new SGJoin(prevHook,
                                    SGSide(SGPosition(newSeq->getID(), 0),
                                           false));
        _sg->addJoin(inJoin);
      }
      
      if (j < dnaLength - 1)
      {
        SGSide tgtSide(sgPositions[j+1], false);
        SGJoin* outJoin = new SGJoin(
          SGSide(SGPosition(newSeq->getID(), newSeq->getLength() -1), true),
          tgtSide);
        _sg->addJoin(outJoin);
      }

      prevHook.setBase(SGPosition(newSeq->getID(), newSeq->getLength() - 1));
      prevHook.setForward(true);

      for (sg_int_t k = i; k <= j; ++k)
      {
        sgPositions[k].setSeqID(newSeq->getID());
        sgPositions[k].setPos(k - i);

        cout << "addSNP " << SGPosition(tgtPos.getSeqID(), tgtPos.getPos() + k)
             << " , " <<  dnaString[dnaOffset + k] << " , "
             << sgPositions[k] << endl;
        
        addSNP(SGPosition(tgtPos.getSeqID(), tgtPos.getPos() + k),
               dnaString[dnaOffset + k],
               sgPositions[k]);
      }
      
      i = j;
    }
    else
    {
      prevHook.setBase(sgPositions[i]);
      prevHook.setForward(true);
    }
  }

  // now the side graph is updated, and sgPositions 
  // finally we update the srcLookup for every snp, to make sure that
  // we always map to the right base (ie when using the lookup structure
  // to compute paths for the input

  for (sg_int_t i = 0; i < dnaLength; ++i)
  {
    sg_int_t j = i + 1;
    while (j < dnaLength && sgPositions[j].getSeqID() ==
           sgPositions[j-1].getSeqID() &&
           sgPositions[j].getPos() ==
           sgPositions[j-1].getPos() + 1)
    {
      ++j;
    }
    SGPosition pos = srcPos;
    pos.setPos(srcPos.getPos() + i);
    srcLookup->addInterval(pos, sgPositions[i], j - i, false);
  }

  // return hooks at either end to be joined by calling code to rest of
  // graph
  pair<SGSide, SGSide> outHooks;
  outHooks.first.setBase(sgPositions[0]);
  outHooks.first.setForward(false);
  outHooks.second.setBase(sgPositions[sgPositions.size() - 1]);
  outHooks.second.setForward(true);
  
  return outHooks;
}

const SGPosition& SNPHandler::findSNP(const SGPosition& pos, char nuc)
{
  if (_caseSens == false)
  {
    nuc = toupper(nuc);
  }
  
  if (pos != _cachePos)
  {
    pair<SNPMap::iterator, bool> ins = _snpMap.insert(
      pair<SGPosition, SNPList*>(pos, NULL));
    _cacheIt = ins.first;
    _cachePos = pos;
  }

  SNPList* snpList = _cacheIt->second;
  if (snpList != NULL)
  {
    for (SNPList::iterator i = snpList->begin(); i != snpList->end(); ++i)
    {
      if (i->_nuc == nuc)
      {
        return i->_pos;
      }
    }
  }
  return SideGraph::NullPos;
}

void SNPHandler::addSNP(const SGPosition& pos, char nuc,
                        const SGPosition& snpPosition)
{
  assert(findSNP(pos, nuc) == SideGraph::NullPos);

  if (_caseSens == false)
  {
    nuc = toupper(nuc);
  }
         
  if (pos != _cachePos)
  {
    pair<SNPMap::iterator, bool> ins = _snpMap.insert(
      pair<SGPosition, SNPList*>(pos, NULL));
    _cacheIt = ins.first;
    _cachePos = pos;
  }
  if (_cacheIt->second == NULL)
  {
    _cacheIt->second = new SNPList();
  }
  
  SNPList* snpList = _cacheIt->second;
  snpList->push_back(SNP(snpPosition, nuc));

  // add new position into the handler
  assert(_snpMap.find(snpPosition) == _snpMap.end());
  _snpMap.insert(pair<SGPosition, SNPList*>(snpPosition, snpList));
}

void SNPHandler::getSNPName(const SGPosition& tgtPos, const string& dnaString,
                            sg_int_t dnaOffset, bool reverseMap, sg_int_t i,
                            string& nameBuf) const
{
  const SGSequence* seq = _sg->getSequence(tgtPos.getSeqID());
  stringstream ss;
  ss << seq->getName() << "_" << (tgtPos.getPos() + dnaOffset + i)
     << "_SNP_" << dnaString[dnaOffset + i];
  nameBuf = ss.str();
}
