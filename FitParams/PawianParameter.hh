//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 Julian Pychy (julian@ep1.rub.de)                       //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

#pragma once

#include <algorithm>
#include <memory>
#include <cassert>
#include <string> 

class PawianParameter {

public:
  
   //constructor for constant Parameter
   PawianParameter(unsigned int num, const std::string & name, double val) : 
      _num(num), _value(val), _error(0.), _const(true), _fix(false), 
      _loLimit(0.), _upLimit(0.), _loLimValid(false), _upLimValid(false),
      _name(name)
   {}
  
   //constructor for standard Parameter
   PawianParameter(unsigned int num, const std::string & name, double val, double err) :
      _num(num), _value(val), _error(err), _const(false), _fix(false), 
      _loLimit(0.), _upLimit(0.), _loLimValid(false), _upLimValid(false),
      _name(name)
   {}
  
   //constructor for limited Parameter
   PawianParameter(unsigned int num, const std::string & name, double val, double err, 
                   double min, double max) : 
      _num(num),_value(val), _error(err), _const(false), _fix(false), 
      _loLimit(min), _upLimit(max), _loLimValid(true), _upLimValid(true), 
      _name(name)    
   {
      assert(min != max);
      if(min > max) {
         _loLimit = max;
         _upLimit = min;
      }
   }

   ~PawianParameter() {}

   PawianParameter(const PawianParameter& par) : 
      _num(par._num), _value(par._value), _error(par._error),
      _const(par._const), _fix(par._fix), _loLimit(par._loLimit), 
      _upLimit(par._upLimit), _loLimValid(par._loLimValid), 
      _upLimValid(par._upLimValid), 
      _name(par._name ) 
   {}

   PawianParameter& operator=(const PawianParameter& par) {
      _num = par._num;
      _name = par._name;
      _value = par._value;
      _error = par._error;
      _const = par._const;
      _fix = par._fix;
      _loLimit = par._loLimit; 
      _upLimit = par._upLimit;
      _loLimValid = par._loLimValid; 
      _upLimValid = par._upLimValid;
      return *this;
   }

   //access methods
   unsigned int Number() const {return _num;}
   // new API returning a string 
   const std::string & GetName() const { return _name; }
   // return const char * for mantaining backward compatibility
   const char * Name() const {return _name.c_str();}
 
   double Value() const {return _value;}
   double Error() const {return _error;}

   //interaction
   void SetValue(double val) {_value = val;}
   void SetError(double err) {_error = err;}
   void SetLimits(double low, double up) {
      assert(low != up);
      _loLimit = low; 
      _upLimit = up;
      _loLimValid = true; 
      _upLimValid = true;
      if(low > up) {
         _loLimit = up; 
         _upLimit = low;
      }
   }

   void SetUpperLimit(double up) {
      _loLimit = 0.; 
      _upLimit = up;
      _loLimValid = false; 
      _upLimValid = true;
   }

   void SetLowerLimit(double low) {
      _loLimit = low; 
      _upLimit = 0.;
      _loLimValid = true; 
      _upLimValid = false;
   }

   void RemoveLimits() {
      _loLimit = 0.; 
      _upLimit = 0.;
      _loLimValid = false; 
      _upLimValid = false;
   }

   void Fix() {_fix = true;}
   void Release() {_fix = false;}
  
   //state of Parameter (fixed/const/limited)
   bool IsConst() const {return _const;}
   bool IsFixed() const {return _fix;}

   bool HasLimits() const {return _loLimValid || _upLimValid; }
   bool HasLowerLimit() const {return _loLimValid; }
   bool HasUpperLimit() const {return _upLimValid; }
   double LowerLimit() const {return _loLimit;}
   double UpperLimit() const {return _upLimit;}

private:

   unsigned int _num;
   double _value;
   double _error;
   bool _const;
   bool _fix;
   double _loLimit; 
   double _upLimit;
   bool _loLimValid; 
   bool _upLimValid;
   std::string _name;

private:

};


