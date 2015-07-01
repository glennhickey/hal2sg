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
using namespace hal;

SNPHandler::SNPHandler(SideGraph* sideGraph, bool caseSensitive,
                       bool onlySequenceNames)
  :  _caseSens(caseSensitive), _sg(sideGraph), _snpCount(0),
     _onlySequenceNames(onlySequenceNames)
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

pair<SGSide, SGSide> SNPHandler::createSNP(const string& srcDNA,
                                           const string& tgtDNA,
                                           size_t dnaOffset,
                                           size_t dnaLength,
                                           const Sequence* halSrcSequence,
                                           const SGPosition& srcPos,
                                           const SGPosition& tgtPos,
                                           bool reverseMap,
                                           SGLookup* srcLookup,
                                           SGLookBack* seqMapBack)
{
  bool say = false;
  if (say)
  {
  cout << "\ncreateSNP " << srcDNA << " " << dnaOffset << "," << dnaLength
       << ",r=" << reverseMap  << endl;
  cout << "          " << tgtDNA << endl;
  cout << "srange " << srcDNA.substr(dnaOffset, dnaLength) << endl;
  cout << "trange " << tgtDNA.substr(dnaOffset, dnaLength) << endl;
  }
  vector<SGPosition> sgPositions(dnaLength);
  // note : should use hints / cache to speed up consecutive bases
  // but start simple for baseline tests

  // first we use the snp structure to find out if any of the SNPs we
  // want to add already exist.  If they do, we update sgPositions
  SGPosition tgt(tgtPos);
  for (sg_int_t i = 0; i < dnaLength; ++i)
  {
    // tgtDNA is reversed already, but "tgt" coordinates are forward
    // and need to move in right-to-left direction. 
    sg_int_t tgtDelta = !reverseMap ? i : -i;
    tgt.setPos(tgtPos.getPos() + tgtDelta);

    char srcVal = srcDNA[i + dnaOffset];
    char tgtVal = tgtDNA[i + dnaOffset];       
    assert(isSub(srcVal,tgtVal));
    // tgtVal in its forward orientation (this is what will be in the
    // snphandler structure)
    char tgtValForward = !reverseMap ? tgtVal : reverseComplement(tgtVal);
    
    // add a baseline snp for the target (forward dir) if not there already. 
    // create / update any other structures. 
    if (findSNP(tgt, tgtValForward) == SideGraph::NullPos)
    {
      if (say)
      {
       cout << "addtgtSNP " << SGPosition(tgtPos.getSeqID(), tgt.getPos())
            << " , " <<  srcVal << "->" << tgtValForward << " , "
              << SGPosition(tgtPos.getSeqID(), tgt.getPos()) << endl;
      }
      assert(srcVal != tgtVal);
      addSNP(SGPosition(tgtPos.getSeqID(), tgt.getPos()),
             tgtValForward, 
             SGPosition(tgtPos.getSeqID(), tgt.getPos()));
    }

    if (say)
    {
      cout << "findSNP " << tgt << ","
           << "srcVal=" << srcVal << "tgtVal=" << tgtVal
           << " (srcpos="
           << srcPos.getPos() + i << ") " 
          << " = " << findSNP(tgt, srcVal)
          << " i=" << i << ", tgtPos=" << tgtPos << " rm " << reverseMap << endl;
    }
    sgPositions[i] = findSNP(tgt, srcVal);
  }

  // now we have a position for all existing snps in the graph.  any
  // positions that are "null" means we need to add new sequences to
  // cover them.  We create the sequences here and add them
  // to the side graph.  We also update sgPositions with these new
  // coordinates.
  set<sg_int_t> newSnpSeqs;
  SGSide prevHook;
  string nameBuf;
  
  for (sg_int_t i = 0; i < dnaLength; ++i)
  {
    if (sgPositions[i] == SideGraph::NullPos)
    {
      sg_int_t j = i + 1;
      for (; j < dnaLength && sgPositions[j] == SideGraph::NullPos; ++j);
      --j;
      getSNPName(halSrcSequence, srcPos, i, j - i + 1, reverseMap, nameBuf);
      const SGSequence* newSeq;
      newSeq = _sg->addSequence(new SGSequence(-1, j - i + 1, nameBuf));
      _snpCount += newSeq->getLength();
      newSnpSeqs.insert(newSeq->getID());
            
      prevHook.setBase(SGPosition(newSeq->getID(), newSeq->getLength() - 1));
      prevHook.setForward(true);

      for (sg_int_t k = i; k <= j; ++k)
      {
        sgPositions[k].setSeqID(newSeq->getID());
        sgPositions[k].setPos(k - i);
        sg_int_t tgtCoord = tgtPos.getPos() + (!reverseMap ? k : -k);
        if (say)
        {
         cout << "addSNP " << SGPosition(tgtPos.getSeqID(), tgtCoord)
              << " , " <<  srcDNA[dnaOffset + k] << " , "
              << sgPositions[k] << endl;
        }
        addSNP(SGPosition(tgtPos.getSeqID(), tgtCoord),
               srcDNA[dnaOffset + k],
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

    if (say)
    {
    cout << "add snp map " << pos << " --> " << sgPositions[i]
         << " len=" << (j-1) << " rev=" << false << endl;
    }
    if (seqMapBack != NULL &&
        newSnpSeqs.find(sgPositions[i].getSeqID()) != newSnpSeqs.end())
    {
      // keep record of where it came from (ie to trace back from the side
      // graph to hal (only made optional to let unit tests skip this step
      // if needed)
      seqMapBack->addInterval(
        SGPosition(sgPositions[i].getSeqID(), 0),
        halSrcSequence,
        pos.getPos(),
        _sg->getSequence(sgPositions[i].getSeqID())->getLength(),
        false);

      if (say)
      {
      cout << "add snp back " << SGPosition(sgPositions[i].getSeqID(), 0)
           << " --> " << halSrcSequence->getFullName() << ", "
           << pos.getPos() << " len="
           << _sg->getSequence(sgPositions[i].getSeqID())->getLength()
           << " rev=false" << endl;
      }
    }
    i = j - 1;
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
  if (pos != snpPosition)
  {
    assert(_snpMap.find(snpPosition) == _snpMap.end());
    _snpMap.insert(pair<SGPosition, SNPList*>(snpPosition, snpList));
  }
}

void SNPHandler::getSNPName(const Sequence* halSrcSequence,
                            const SGPosition& srcPos,
                            sg_int_t offset, sg_int_t length,
                            bool reverseMap, string& outName) const
{
  stringstream ss;
  if (halSrcSequence != NULL)
  {
    // unit tests sometimes dont bother with a hal sequence so we
    // let it be optional here. 
    ss << (_onlySequenceNames ? halSrcSequence->getName() :
           halSrcSequence->getFullName());
  }
  ss << "_" << (srcPos.getPos() + offset) << "_" << length;
  outName = ss.str();
}
