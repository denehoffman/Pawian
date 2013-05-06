// TensorIndexOrder class definition file -*- C++ -*-
/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
// Author: Mike Williams
// Author: Matt Shepherd
#ifndef _TensorIndexOrder_h
#define _TensorIndexOrder_h
//_____________________________________________________________________________

#include <cassert>
#include <iostream>

using namespace std;
//_____________________________________________________________________________
/** @file TensorIndexOrder.h
 *  @brief TensorIndexOrder class definition file.
 */
//_____________________________________________________________________________

#include <iostream>
class TensorIndexOrder {

private:

  unsigned int _rank; ///< Tensor rank this deals with
  int _index; ///< Index to Tensor vector of elements

public:

  // create/copy/destroy
  TensorIndexOrder() : _rank(0),_index(0){/** Default Constructor */}

  /// Constructor (connect with rank @a rank tensor)
  TensorIndexOrder(unsigned int __rank) : _index(0){		
    assert( __rank < 8 );
    _rank = __rank;
  }
    
  virtual ~TensorIndexOrder() {/** Destructor */}

  // operators:

  /// Returns the tensor index value for @a level
  inline int operator[](unsigned int __level) const {
    return ((_index >> (__level << 2 ) ) & 0xF);
  }

  // permutation functions:

  /** Obtains the next valid permuation of the current stored index.
   *
   * This function is NOT to be used when the TensorIndexOrder object is being
   * used to get elements from a tensor (the normal use). 
   * This function is meant to be used on a TensorIndexOrder object declared for 
   * the sole purpose of performing permutations on the indicies of a Tensor.
   */
  TensorIndexOrder& Permute();

  /** Checks to see if the current index permutation is valid.
   *
   * Valid means all indicies are between 0 and Rank-1. The function doesn't 
   * actually check if any 2 indicies are equal, this is done in 
   * TensorIndexOrder::Permutation.
   */
  bool PermIsValid() const {
    for(int i = 0; i < this->Size(); i++){
      if(((*this)[i] < 0)||((*this)[i] >= this->Size())) return false;
    }
    return true;
  }

  // functions:

  /// Returns the number of entries (coresponds to tensor rank)
  inline int Size() const {return _rank;}

  /// Set the index value
  inline void SetIndex(unsigned int __i,unsigned int __value){
     //    assert(__value < 16 && __i < _rank );
    _index = (_index & ~(0xF << (__i << 2))) + (__value << (__i << 2));
  }
};
//_____________________________________________________________________________

#endif /* _TensorIndexOrder_H */
