/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SGSEQUENCE_H
#define _SGSEQUENCE_H


#include "sgcommon.h"

/**
 * Named sequence.  We will rely on separate map to/from HAL sequences
 * in order to get DNA, I think...  IDs will be unique numbers from
 * 0,1,2... since no plans to do anything but add, at this point.  
 */
class SGSequence
{
public:

   SGSequence();
   SGSequence(sg_seqid_t id, sg_int_t length, const std::string& name);
   ~SGSequence();

   sg_seqid_t getID() const;
   sg_int_t getLength() const;
   const std::string& getName() const;

   void set(sg_seqid_t id, sg_int_t length, const std::string& name);
   void setID(sg_seqid_t id);
   void setLength(sg_int_t length);
   void setName(const std::string& name);

protected:

   sg_seqid_t _id;
   sg_int_t _length;
   std::string _name;
};

std::ostream& operator<<(std::ostream& os, const SGSequence& s);

inline SGSequence::SGSequence() : _id(-1), _length(-1)
{
}

inline SGSequence::SGSequence(sg_seqid_t id, sg_int_t length,
                              const std::string& name) :
  _id(id), _length(length), _name(name)
{
}

inline SGSequence::~SGSequence()
{
}

inline sg_seqid_t SGSequence::getID() const
{
  return _id;
}

inline sg_int_t SGSequence::getLength() const
{
  return _length;
}

inline const std::string& SGSequence::getName() const
{
  return _name;
}

inline void SGSequence::set(sg_seqid_t id, sg_int_t length,
                            const std::string& name)
{
  _id = id;
  _length = length;
  _name = name;
}

inline void SGSequence::setID(sg_seqid_t id)
{
  _id = id;
}

inline void SGSequence::setLength(sg_int_t length)
{
  _length = length;
}

inline void SGSequence::setName(const std::string& name)
{
  _name = name;
}

inline std::ostream& operator<<(std::ostream& os, const SGSequence& s)
{
  return os << "s(" << s.getID() << "," << s.getLength() << ","
            << s.getName() << ")";
}

#endif
