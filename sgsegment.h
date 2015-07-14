/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGSEGMENT_H
#define _SGSEGMENT_H

#include "sgside.h"

/**
* The length of this Segment's sequence. If start is on the forward strand, 
* the Segment contains the range 
* [start.base.position, start.base.position + length). If start is on the 
* reverse strand, the Segment contains the range 
* (start.base.position - length, start.base.position]. This is equivalent to
* starting from the side indicated by start, and traversing through that base
*  out to the specified length.
*
* A Segment may have zero length (for example, when it is being used to 
* specify a Path consisting only of a Join.
*/ 
class SGSegment
{
public:

   SGSegment();
   SGSegment(const SGSide& side, sg_int_t length);
   SGSegment(const SGSide& start, const SGSide& end);
   ~SGSegment();

   const SGSide& getSide() const;
   sg_int_t getLength() const;

   void set(const SGSide& side, sg_int_t length);
   void setSide(const SGSide& side);
   void setLength(sg_int_t length);
   /** flip the direction of the segment (but it will still cover same
    * positions)
    */
   void flip();

   /** get the minimum value endpoint
    */
   SGPosition getMinPos() const;

   /** get the maximum value endpoint
    */
   SGPosition getMaxPos() const;

   /** get the side we'd join to if we're going into the segment
    */
   SGSide getInSide() const;

   /** get the side we'd join from if we're going out of the segment
    */
   SGSide getOutSide() const;
      
   bool operator<(const SGSegment& s2) const;
   bool operator<=(const SGSegment& s2) const;
   bool operator==(const SGSegment& s2) const;
   bool operator!=(const SGSegment& s2) const;
   
protected:
   
   SGSide _side;
   sg_int_t _length;
};
std::ostream& operator<<(std::ostream& os, const SGSegment& s);

inline SGSegment::SGSegment() : _length(-1)
{
}

inline SGSegment::SGSegment(const SGSide& side, sg_int_t length) :
  _side(side), _length(length)
{
}

inline SGSegment::SGSegment(const SGSide& start, const SGSide& end) :
  _side(start)
{
    _length = start.lengthTo(end);
    
  // Segment definition is always [x,y) so we have to adjust our start
  // accordingly
  const SGPosition& startPos = start.getBase();
  if (!start.getForward() && end >= start)
  {
    _side.setBase(SGPosition(startPos.getSeqID(), startPos.getPos() + 1));
    _side.setForward(true);
  }
  else if (start.getForward() && end <= start)
  {
    _side.setBase(SGPosition(startPos.getSeqID(), startPos.getPos() - 1));
    _side.setForward(false);
  }
}

inline SGSegment::~SGSegment()
{
}

inline const SGSide& SGSegment::getSide() const
{
  return _side;
}

inline sg_int_t SGSegment::getLength() const
{
  return _length;
}

inline void SGSegment::set(const SGSide& side, sg_int_t length)
{
  _side = side;
  _length = length;
}

inline void SGSegment::setSide(const SGSide& side)
{
  _side = side;
}

inline void SGSegment::setLength(sg_int_t length)
{
  _length = length;
}

inline void SGSegment::flip()
{
  if (_length > 0)
  {
    sg_int_t delta = _side.getForward() ? _length - 1 : -_length + 1;
    _side.setBase(SGPosition(_side.getBase().getSeqID(),
                             _side.getBase().getPos() + delta));
  }
  _side.setForward(!_side.getForward());
}

inline SGPosition SGSegment::getMinPos() const
{
  return _side.getForward() == true ? _side.getBase() :
     SGPosition(_side.getBase().getSeqID(),
                _side.getBase().getPos() - _length + 1);
}

inline SGPosition SGSegment::getMaxPos() const
{
  return _side.getForward() == false ? _side.getBase() :
     SGPosition(_side.getBase().getSeqID(),
                _side.getBase().getPos() + _length - 1);
}

inline SGSide SGSegment::getInSide() const
{
  // if we're forward, we join on the (left) forward
  return _side;
}

inline SGSide SGSegment::getOutSide() const
{
  SGSide side = _side;
  if (_length > 0)
  {
    sg_int_t delta = _side.getForward() ? _length - 1 : -_length + 1;
    side.setBase(SGPosition(_side.getBase().getSeqID(),
                            _side.getBase().getPos() + delta));
  }
  // if we're forward, want to join to reverse
  side.setForward(!_side.getForward());
  return side;
}

inline bool SGSegment::operator<(const SGSegment& s2) const
{
  return _side < s2._side || (_side == s2._side &&
                              _length < s2._length);
}

inline bool SGSegment::operator<=(const SGSegment& s2) const
{
  return *this < s2 || *this == s2;
}

inline bool SGSegment::operator==(const SGSegment& s2) const
{
  return _side == s2._side && _length == s2._length;
}

inline bool SGSegment::operator!=(const SGSegment& s2) const
{
  return !(*this == s2);
}

inline std::ostream& operator<<(std::ostream& os, const SGSegment& s)
{
  return os << "seg(" << s.getSide() << "," << s.getLength() << ")";
}

#endif
