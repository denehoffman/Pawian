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
/** @class TensorIndexOrder
 *  @author Mike Williams
 *  @author Matt Shepherd
 *
 *  @brief A class that converts \f$(\mu,\nu,...)\f$ to an index for Tensor.
 *
 * This class is a utility class used by Tensor. Its current implementation is
 * limited to tensors of rank 15 or less (well above current memory 
 * capacities).
 * The elements of Tensor<type> T(rank) can be accessed via           
 * TensorIndexOrder index(rank) by simply calling T(index). 
 *
 * <b> Example Usage </b>
 *
 * <!--
 *      TensorIndexOrder index(this->Rank());              
 *      while(index.IsValid()){                       
 *        (*this)(index) "some operations ...";       
 *        ++index;                                                      
 *      }                                                              
 * -->
 * \include TensorIndexOrder.ex
 *
 * Note: Indicies are t = 0, x = 1, y = 2 and z = 3 
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

  /// Returns the current index
  inline int operator()() const {return _index;}

  /** Prefix step up incrementer.
   *
   * If TensorIndexOrder index(rank) is defined, then ++index steps up index
   * to the next valid tensor index.
   *
   * Example: If TensorIndexOrder index(3) = (0,2,1), then ++index = (0,2,2)). 
   *
   * note: Postfix incrementer, index++, is also defined.
   */
  TensorIndexOrder& operator++(){ _index++; return (*this); }

  /// Postfix step up incrementer (see operator++()).
  inline TensorIndexOrder& operator++(int){
    return (++(*this));
  }

  /** Prefix step down incrementer.
   *
   * If TensorIndexOrder index(rank) is defined, then --index steps down index
   * to the next valid tensor index.
   *
   * Example: For TensorIndexOrder index(3) = (0,2,1), --index = (0,2,0). 
   *
   * note: Postfix incrementer, index--, is also defined.
   */
  TensorIndexOrder& operator--(){ _index--; return (*this);}

  /// Postfix step down incrementer (see operator--()).
  inline TensorIndexOrder& operator--(int){
    return (--(*this));
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

  /// Adjust the rank
  inline void Resize(unsigned int __rank) {     
     std::cout << "rcreated tio " << __rank << std::endl;
    assert(__rank < 8 ); 
    _rank = __rank; 
  }

  /// Set the index value
  inline void SetIndex(unsigned int __i,unsigned int __value){
     //    assert(__value < 16 && __i < _rank );
    _index = (_index & ~(0xF << (__i << 2))) + (__value << (__i << 2));
  }

  /// Checks to see if the current tensor index stored is a valid index.
  inline bool IsValid() const {
    if( ( _index >> ( _rank << 2 ) ) != 0 ) return false;
    return true;
  }

  /// Resets all indicies in this object to 0.
  inline void Reset() {_index = 0;}

  /** Print to screen the current set of indicies.
   *  Format is: \f$ (\mu,\nu,\rho,...,\pi) \f$ 
   */
  void Print(ostream &__os = cout);
};
//_____________________________________________________________________________

#endif /* _TensorIndexOrder_H */
