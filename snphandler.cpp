/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <limits>
#include <cctype>
#include <cstdio>

#include "snphandler.h"

using namespace std;

SNPHandler::SNPHandler(bool caseSensitive) : _caseSens(caseSensitive)
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

const SGPosition& SNPHandler::findSNP(const SGPosition& pos, char nuc)
{
  if (_caseSens == false)
  {
    nuc = std::toupper(nuc);
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
    nuc = std::toupper(nuc);
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
  _snpMap.insert(pair<SGPosition, SNPList*>(snpPosition, snpList));
}
