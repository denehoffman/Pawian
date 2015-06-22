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

// PwaCovMatrix class definition file. -*- C++ -*-
// Copyright 2015 Bertram Kopf

#pragma once

#include <iostream> 
#include <vector>
#include <memory>
#include <map>
#include "Minuit2/MnUserParameters.h"

using namespace ROOT::Minuit2;

//class MnUserParameters;

class AbsPawianParameters {

public:

  AbsPawianParameters();

  ~AbsPawianParameters();

  virtual AbsPawianParameters& operator=(AbsPawianParameters& par)=0;
  virtual AbsPawianParameters* Clone()=0;

  virtual std::string type()=0;

  virtual std::shared_ptr<MnUserParameters> mnUserParametersPtr();


  virtual unsigned int VariableParameters() const=0;

  virtual const std::vector<std::string> ParamNames()=0;

   // access to parameters and errors in column-wise representation 
   virtual std::vector<double> Params() const=0;
   virtual std::vector<double> Errors() const=0;

   // access to single Parameter

   // Add free Parameter Name, Value, Error
   virtual bool Add(const std::string &, double, double)=0;
   // Add limited Parameter Name, Value, Lower bound, Upper bound
   virtual bool Add(const std::string &, double, double, double, double)=0;
   // Add const Parameter Name, value
   virtual bool Add(const std::string &, double)=0;

   /// interaction via external number of Parameter
   virtual void Fix(unsigned int)=0;
   virtual void Release(unsigned int)=0;
   virtual void SetValue(unsigned int, double)=0;
   virtual void SetError(unsigned int, double)=0;
   virtual void SetLimits(unsigned int, double, double)=0;
   virtual void SetUpperLimit(unsigned int, double)=0;
   virtual void SetLowerLimit(unsigned int, double)=0;
   virtual void RemoveLimits(unsigned int)=0;

   virtual double Value(unsigned int) const=0;
   virtual double Error(unsigned int) const=0;
  
   /// interaction via Name of Parameter
   virtual void Fix(const std::string &)=0;
   virtual void Release(const std::string &)=0;
   virtual void SetValue(const std::string &, double)=0;
   virtual void SetError(const std::string &, double)=0;
   virtual void SetLimits(const std::string &, double, double)=0;
   virtual void SetUpperLimit(const std::string &, double)=0;
   virtual void SetLowerLimit(const std::string &, double)=0;
   virtual void RemoveLimits(const std::string &)=0;

   virtual double Value(const std::string &) const=0;
   virtual double Error(const std::string &) const=0;

  /// interaction via external number of Parameter
  virtual bool IsFixed(unsigned int)=0;
  virtual bool HasLimits(unsigned int)=0;

  virtual double UpperLimit(unsigned int)=0;
  virtual double LowerLimit(unsigned int)=0; 
 
   //convert Name into external number of Parameter
   virtual unsigned int Index(const std::string &) const=0;
   //convert external number into Name of Parameter
   virtual const std::string & GetName(unsigned int) const=0;
   // mantain interface with const char * for backward compatibility
   virtual const char* Name(unsigned int) const=0;

  virtual void print(std::ostream& os, bool extended=false);
  virtual void SetAllValues(const std::vector<double>& values);
  virtual std::shared_ptr<AbsPawianParameters> paramsWithSameOrder(std::shared_ptr<AbsPawianParameters> toSort);
  virtual void SetDependencies(const std::string &, std::vector<unsigned int>& depVec);
  virtual void SetDependencies(unsigned int, std::vector<unsigned int>& depVec);

protected:
  std::map<unsigned int, std::vector<unsigned int> > _depMap;

private:
};

