/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGPOSITION_H
#define _SGPOSITION_H


#include "sgsequence.h"

/**
 *  Position in the sidegraph
 */
class SGPosition
{
public:
   
   SGPosition();
   SGPosition(sg_seqid_t seqID, sg_int_t pos);
   ~SGPosition();

   sg_seqid_t getSeqID() const;
   sg_int_t getPos() const;

   void set(sg_seqid_t seqID, sg_int_t pos);
   void setSeqID(sg_seqid_t seqID);
   void setPos(sg_int_t pos);

   bool operator<(const SGPosition& p2) const;
   bool operator==(const SGPosition& p2) const;
   bool operator!=(const SGPosition& p2) const;
   
protected:     

   sg_seqid_t _seqid;
   sg_int_t _pos;
};

std::ostream& operator<<(std::ostream& os, const SGPosition& p);


inline SGPosition::SGPosition() : _seqid(-1), _pos(-1)
{
}

inline SGPosition::SGPosition(sg_seqid_t seqID, sg_int_t pos) :
  _seqid(seqID), _pos(pos)
{
}

inline SGPosition::~SGPosition()
{
}

inline sg_seqid_t SGPosition::getSeqID() const
{
  return _seqid;
}

inline sg_int_t SGPosition::getPos() const
{
  return _pos;
}

inline void SGPosition::set(sg_seqid_t seqID, sg_int_t pos)
{
  _seqid = seqID;
  _pos = pos;
}

inline void SGPosition::setSeqID(sg_seqid_t seqID)
{
  _seqid = seqID;
}

inline void SGPosition::setPos(sg_int_t pos)
{
  _pos = pos;
}

inline bool SGPosition::operator<(const SGPosition& p2) const
{
  return _seqid < p2._seqid || (_seqid == p2._seqid && _pos < p2._pos);
}

inline bool SGPosition::operator==(const SGPosition& p2) const
{
  return _seqid == p2._seqid && _pos == p2._pos;
}

inline bool SGPosition::operator!=(const SGPosition& p2) const
{
  return _seqid != p2._seqid || _pos != p2._pos;
}

inline std::ostream& operator<<(std::ostream& os,
                                const SGPosition& p)
{
  return os << "p(" << p.getSeqID() << "," << p.getPos() << ")";
}

#endif
