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
  for (SNPMap::iterator i = _snpMap.begin(); i != _snpMap.end(); ++i)
  {
    delete i->second;
  }
}

const SGPosition& SNPHandler::findSNP(const SGPosition& pos, char nuc)
{
  if (pos != _cachePos)
  {
    pair<SNPMap::iterator, bool> ins = _snpMap.insert(
      pair<SGPosition, SNPList*>(pos, NULL));
    _cacheIt = ins.first;
    _cachePos = pos;
  }

  if (_caseSens == false)
  {
    nuc = std::toupper(nuc);
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
  if (pos != _cachePos)
  {
    pair<SNPMap::iterator, bool> ins = _snpMap.insert(
      pair<SGPosition, SNPList*>(pos, NULL));
    _cacheIt = ins.first;
    _cachePos = pos;
    if (_cacheIt->second == NULL)
    {
      _cacheIt->second = new SNPList();
    }
  }
  SNPList* snpList = _cacheIt->second;

  if (_caseSens == false)
  {
    nuc = std::toupper(nuc);
  }

#ifndef _NDEBUG
  for (SNPList::iterator i = snpList->begin(); i != snpList->end(); ++i)
  {
    if (i->_nuc == nuc)
    {
      assert(i->_nuc != nuc);
    }
  }
#endif

  snpList->push_back(SNP(snpPosition, nuc));
  
}
