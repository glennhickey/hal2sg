/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGSIDE_H
#define _SGSIDE_H

#include "sgposition.h"

/**
 * Position of join endpoint
 */ 
class SGSide
{
public:

   SGSide();
   SGSide(const SGPosition& base, bool forward);
   ~SGSide();

   const SGPosition& getBase() const;
   bool getForward() const;

   void set(const SGPosition& base, bool forward);

   bool operator<(const SGSide& s2) const;
   bool operator==(const SGSide& s2) const;
   
protected:
   
   SGPosition _base;
   bool _forward;
};
std::ostream& operator<<(std::ostream& os, const SGSide& s);

inline SGSide::SGSide() : _forward(true)
{
}

inline SGSide::SGSide(const SGPosition& base, bool forward) :
  _base(base), _forward(forward)
{
}

inline SGSide::~SGSide()
{
}

inline const SGPosition& SGSide::getBase() const
{
  return _base;
}

inline bool SGSide::getForward() const
{
  return _forward;
}

inline void SGSide::set(const SGPosition& base, bool forward)
{
  _base = base;
  _forward = forward;
}

inline bool SGSide::operator<(const SGSide& s2) const
{
  return _base < s2._base || (_base == s2._base &&
                              _forward < s2._forward);
}

inline bool SGSide::operator==(const SGSide& s2) const
{
  return _base == s2._base && _forward == s2._forward;
}

inline std::ostream& operator<<(std::ostream& os, const SGSide& s)
{
  return os << "s(" << s.getBase() << "," << s.getForward() << ")";
}

#endif
