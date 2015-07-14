/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGJOIN_H
#define _SGJOIN_H

#include "sgside.h"

/**
 * Join from one side to another
 */
class SGJoin
{
public:
   
   SGJoin();
   SGJoin(const SGSide& side1, const SGSide& side2);
   ~SGJoin();
   
   const SGSide& getSide1() const;
   const SGSide& getSide2() const;
   
   void set(const SGSide& side1, const SGSide& side2);
   void setSide1(const SGSide& side1);
   void setSide2(const SGSide& side2);
   /** swap side1 with side2 */
   void swap();
   /** joins two consecutive sides */
   bool isTrivial() const;
   
   bool operator<(const SGJoin& j2) const;
   bool operator==(const SGJoin& j2) const;

protected:
   
   SGSide _side1;
   SGSide _side2;
};

std::ostream& operator<<(std::ostream& os, const SGJoin& j);

struct SGJoinPtrLess
{
   bool operator()(const SGJoin* j1, const SGJoin* j2) const {
     return *j1 < *j2;
   }
};


inline SGJoin::SGJoin()
{
}

inline SGJoin::SGJoin(const SGSide& side1, const SGSide& side2) :
  _side1(side1), _side2(side2)
{
}

inline SGJoin::~SGJoin()
{
}

inline const SGSide& SGJoin::getSide1() const
{
  return _side1;
}

inline const SGSide& SGJoin::getSide2() const
{
  return _side2;
}

inline void SGJoin::set(const SGSide& side1, const SGSide& side2)
{
  _side1 = side1;
  _side2 = side2;
}

inline void SGJoin::setSide1(const SGSide& side1)
{
  _side1 = side1;
}

inline void SGJoin::setSide2(const SGSide& side2)
{
  _side2 = side2;
}

inline void SGJoin::swap()
{
  std::swap(_side1, _side2);
}

inline bool SGJoin::isTrivial() const
{
  return _side1.lengthTo(_side2) == 0 &&
     _side1.getForward() != _side2.getForward();
}

inline bool SGJoin::operator<(const SGJoin& j2) const
{
  return _side1 < j2._side1 || (_side1 == j2._side1 && _side2 < j2._side2);
}

inline bool SGJoin::operator==(const SGJoin& j2) const
{
  return _side1 == j2._side1 && _side2 == j2._side2;
}

inline std::ostream& operator<<(std::ostream& os, const SGJoin& j) {
  return os << "j(" << j.getSide1() << "," << j.getSide2() << ")";
}

#endif
