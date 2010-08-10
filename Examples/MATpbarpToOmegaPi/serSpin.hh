// Spin class definition file -*- C++ -*-
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
#ifndef _serSpin_H
#define _serSpin_H

#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/text_iarchive.hpp>

#include "qft++/relativistic-quantum-mechanics/Spin.hh"

using namespace std;
//_____________________________________________________________________________
/** @file Spin.h
 *  @brief Spin class definition file
 */
//_____________________________________________________________________________
/** @class Spin
 *  @author Mike Williams
 *
 *  @brief A simple class for handling integral and half-integral numbers.
 *
 * The Spin class provides a converion opeator to double, thus it can be used
 * anywhere a double can be used. Spin objects can only take on integral or
 * half-integral values. Any attempt to set a spin object to some other value
 * results in the program being aborted.
 *
 */
//_____________________________________________________________________________

class serSpin {

  friend class boost::serialization::access;

  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) {
    using boost::serialization::make_nvp;

    ar & BOOST_SERIALIZATION_NVP(_numer);
    ar & BOOST_SERIALIZATION_NVP(_denom);

    /* Add your own class-variables here in the following way:
      ar & BOOST_SERIALIZATION_NVP(myVar);
      or
      ar & make_nvp("myVar", myVar); // The latter form can be necessary when dealing with templates
    */
  }

private:
  // data members:
  int _numer; ///< numerator
  int _denom; ///< denominator

  // private functions:

  /// Copy @a spin data members
  void _Copy(const serSpin &__spin){
    _numer = __spin._numer;
    _denom = __spin._denom;
  }

  /// Build a rational spin from a double
  void _BuildFromDouble(double __x){
    int x_int = (int)__x;
    if(abs(__x - x_int) < 1.e-5) { // integer spin
      _numer = x_int;
      _denom = 1;
    }
    else{
      double x = 2.*__x;
      x_int = (int)x;
      if(abs(x - x_int) < 1.e-5) { // half-integer spin
	_numer = x_int;
	_denom = 2;
      }
      else{ // illegal assignment attempt
	cout << "Error! Attempt to assign a Spin object the value " << __x
	     << " which is neither integral nor half-integral." << endl;
	abort();
      }
    }
  }

public:
  
  // create/copy/destroy:

  /// Default Constructor (spin = 0)
  serSpin(){
    _numer = 0;
    _denom = 1;
  }

  /** Default Constructor (spin = @a numerator/denominator)
   * Note: @a denominator must be either 1 or 2 
   */
  serSpin(int __numerator,int __denominator = 1){
    if(__denominator != 1 && __denominator != 2){
      cout << "Error! Attempt to create a spin object with denominator = " 
	   << __denominator << ". Must be either 1 or 2." << endl;
    }
    assert(__denominator == 1 || __denominator == 2);
    _numer = __numerator;
    _denom = __denominator;
  }

  /// Constructor (@a spin must be a valid integral or half-integral spin)
  serSpin(double __spin){
    this->_BuildFromDouble(__spin);
  }

  /// Copy Constructor
  serSpin(const serSpin &__spin){
    this->_Copy(__spin);
  }

  /// Copy Constructor
  serSpin(const Spin &__spin){
    _numer = __spin.Numerator();
    _denom = __spin.Denominator();
  }

  /// Destructor
  virtual ~serSpin(){}

  // Getters:

  /// Returns the numerator of the spin
  inline int Numerator() const {
    return _numer;
  }

  /// Returns the denominator of the spin
  inline int Denominator() const {
    return _denom;
  }

  /// Returns the value (as floating point) of the spin
  inline double Value() const {
    return  _numer/(double)_denom;
  }

  // Setters:

  /// Set the spin
  inline void SetSpin(int __numer,int __denom = 1){
    _numer = __numer;
    _denom = __denom;
  }

  // operators:

  /// Assignment operator
  inline serSpin& operator=(const serSpin &__spin){
    this->_Copy(__spin);
    return *this;
  }

  /// Assignment operator
  inline serSpin& operator=(double __spin){
    this->_BuildFromDouble(__spin);
    return *this;
  }

  /// Conversion operator
  operator double () const {
    return _numer/(double)_denom;
  }

  /// Prefix step up incrementer
  inline serSpin& operator++(){
    _numer += _denom;
    return *this;
  }

  /// Postfix step up incrementer
  inline serSpin& operator++(int){
    return ++(*this);
  }

  /// Prefix step down incrementer
  inline serSpin& operator--(){
    _numer -= _denom;
    return *this;
  }

  /// Postfix step down incrementer
  inline serSpin& operator--(int){
    return --(*this);
  }

  /// Set @a this = @a this + @a spin
  inline serSpin& operator+=(const serSpin &__spin){
    *this = *this + __spin;
    return *this;
  }

  /// Set @a this = @a this - @a spin
  inline serSpin& operator-=(const serSpin &__spin){
    *this = *this - __spin;
    return *this;
  }

  /// Print spin value to @a os
  inline void Print(ostream &__os = cout) const {
    __os << _numer;
    if(_denom != 1) __os << "/" << _denom;
  }
};
//_____________________________________________________________________________

/// Overload of << for Spin class
inline ostream& operator<<(ostream &__os,const serSpin &__spin){
  __spin.Print(__os);
  return __os;
}
//_____________________________________________________________________________

#endif /* _serSpin_H  */
