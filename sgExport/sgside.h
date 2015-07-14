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
 * The forward position is the *left* side of a base. 
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
   void setBase(const SGPosition& base);
   void setForward(bool forward);
   
   /** compute length of semgent that starts with this and 
    * ends with side2 */
   sg_int_t lengthTo(const SGSide& side2) const;
   
   bool operator<(const SGSide& s2) const;
   bool operator<=(const SGSide& s2) const;
   bool operator>=(const SGSide& s2) const;
   bool operator==(const SGSide& s2) const;
   bool operator!=(const SGSide& s2) const;
   
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

inline void SGSide::setBase(const SGPosition& base)
{
  _base = base;
}

inline void SGSide::setForward(bool forward)
{
  _forward = forward;
}

inline sg_int_t SGSide::lengthTo(const SGSide& side2) const
{
  if (this->getBase().getSeqID() != side2.getBase().getSeqID())
  {
    return -1;
  }
  sg_int_t len = -1;
  const SGSide* s1 = this;
  const SGSide* s2 = &side2;
  if (*s2 < *s1)
  {
    std::swap(s1, s2);
  }

  if (s1->getBase() == s2->getBase())
  {
    len = s1->getForward() != s2->getForward() ? 1 : 0;
  }
  else
  {
    len = s2->getBase().getPos() - s1->getBase().getPos() + 1;
    if (s1->getForward() == false)
    {
      --len;
    }
    if (s2->getForward() == true)
    {
      --len;
    }
  }
  return len;
}

inline bool SGSide::operator<(const SGSide& s2) const
{
  return _base < s2._base || (_base == s2._base &&
                              !_forward < !s2._forward);
}

inline bool SGSide::operator<=(const SGSide& s2) const
{
  return *this < s2 || *this == s2;
}

inline bool SGSide::operator>=(const SGSide& s2) const
{
  return !(*this < s2);
}


inline bool SGSide::operator==(const SGSide& s2) const
{
  return _base == s2._base && _forward == s2._forward;
}

inline bool SGSide::operator!=(const SGSide& s2) const
{
  return !(*this == s2);
}

inline std::ostream& operator<<(std::ostream& os, const SGSide& s)
{
  return os << "s(" << s.getBase() << "," << s.getForward() << ")";
}

#endif
