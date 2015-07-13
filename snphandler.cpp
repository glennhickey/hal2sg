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
                                           const SGPosition& sgPos,
                                           bool blockReverseMap,
                                           bool sgReverseMap,
                                           SGLookup* srcLookup,
                                           SGLookBack* seqMapBack)
{
/*   
  cout << "createSNP " << "srcDNA.length=" << srcDNA.length()
       << " dnaOff=" << dnaOffset << " dnaLength=" << dnaLength
       << " srcPos=" << srcPos << " sgPos=" << sgPos
       << " br=" << blockReverseMap << ", sgr=" << sgReverseMap << endl
       << " src = " << srcDNA.substr(dnaOffset, dnaLength) << endl;
       if (blockReverseMap == false)
       cout << " tgt = " << tgtDNA.substr(dnaOffset, dnaLength) << endl;
       else
       cout << " rgt = " << tgtDNA.substr(tgtDNA.length()-1 -(dnaOffset + dnaLength-1), dnaLength) << endl;
*/
   
  assert(srcDNA.length() == tgtDNA.length() && dnaOffset + dnaLength <=
         srcDNA.length());
  // reversal between src and sidegraph (via tgt)
  // recall : blockReverseMap = reversal between src and tgt in block
  //          sgReverseMap = reversal between tgt and sidegraph
  bool tranReverseMap = blockReverseMap != sgReverseMap;

  // run of positions starting from sgPos going forward in side graph
  vector<SGPosition> sgPositions(dnaLength);
  
  // note : should use hints / cache to speed up consecutive bases
  // but start simple for baseline tests

  // first we use the snp structure to find out if any of the SNPs we
  // want to add already exist.  If they do, we update sgPositions
  SGPosition sgCur(sgPos);
  for (sg_int_t i = 0; i < dnaLength; ++i)
  {
    // sgCur is our position in the sidegraph (note that sgPos will be
    // last point in current interface in reversed)
    sg_int_t sgDelta = !tranReverseMap ? i : -i;
    sgCur.setPos(sgPos.getPos() + sgDelta);
    assert(sgCur.getPos() >= 0);

    // srcVal is our corresponding src position
    // (when we reverse mapping, we're actually setting a position
    // in a backwards version that we'll eventually add as new sequence
    // -- very confusing).
    hal_index_t srcIdx =  dnaOffset + i;
    char srcVal = !tranReverseMap ? srcDNA[srcIdx] :
       reverseComplement(srcDNA[srcIdx]);

    // which maps to this position in the target (based on original
    // source position, not reversed srcVal)
    hal_index_t tgtIdx = !blockReverseMap ? dnaOffset + i :
       srcDNA.length() - 1 - dnaOffset - i;
    char tgtVal = tgtDNA[tgtIdx];
    (void)tgtVal;

    // which maps to this position in the sidegraph
    hal_index_t sgIdx = tgtIdx;
    char sgVal = !sgReverseMap ? tgtDNA[sgIdx] :
       reverseComplement(tgtDNA[sgIdx]);
/*
    cout << "  i=" << i <<" srcIDx=" << srcIdx << "," << srcVal
         << " tgtIdx=" << tgtIdx << "," << tgtVal
         << " sgIdx=" << sgIdx << "," << sgVal
         << endl;
*/
  
    // add a baseline snp for (forward) value in the side graph
    if (findSNP(sgCur, sgVal) == SideGraph::NullPos)
    {
      /*
       cout << "addbaseline -> " << sgCur << " = " << sgVal << " -> " << sgCur
         << endl;
      */
      addSNP(sgCur, sgVal, sgCur);
    }

    sgPositions[i] = findSNP(sgCur, srcVal);
    /*
    cout << "sgPositions[" << i << "] = findSnp(" <<sgCur <<","
         <<srcVal <<") =" << sgPositions[i] << endl;
    */
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
      getSNPName(halSrcSequence, srcPos, i, j - i + 1, tranReverseMap, nameBuf);
      const SGSequence* newSeq;
      newSeq = _sg->addSequence(new SGSequence(-1, j - i + 1, nameBuf));
      _snpCount += newSeq->getLength();
      newSnpSeqs.insert(newSeq->getID());

      // note: hooks interface no longer needed -- need to clean everywhere!
      if (!tranReverseMap)
      {
        prevHook.setBase(SGPosition(newSeq->getID(), newSeq->getLength() - 1));
        prevHook.setForward(true);
      }
      else
      {
        prevHook.setBase(SGPosition(newSeq->getID(), 0));
        prevHook.setForward(false);
      }

      for (sg_int_t k = i; k <= j; ++k)
      {
        sg_int_t sgDelta = !tranReverseMap ? k : -j + k - i;
        sgCur.setPos(sgPos.getPos() + sgDelta);

        
        hal_index_t srcIdx = !tranReverseMap ? dnaOffset + k :
           dnaOffset + j - (k - i);
        char srcVal = !tranReverseMap ? srcDNA[srcIdx] :
           reverseComplement(srcDNA[srcIdx]);

        /*
        cout << "k=" << k << " dnaOffset=" << dnaOffset << " i=" << i
             << " j=" << j
             << " dnalength=" << dnaLength << " srcIdx=" << srcIdx
             << " srcVal=" << srcVal << endl;
        */

        // forward map sg to new sequence
        sgPositions[k].setSeqID(newSeq->getID());
        sgPositions[k].setPos(k - i);
        /*
          cout << "addnew -> " << sgCur << " = " << srcVal << " -> "
          << sgPositions[k] << endl;
        */
        addSNP(sgCur, srcVal, sgPositions[k]);
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

    srcLookup->addInterval(pos, sgPositions[i], j - i, tranReverseMap);

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
        tranReverseMap);
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
  assert(pos.getPos() >= 0);
  assert(snpPosition.getPos() >= 0);

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
