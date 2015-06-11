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

#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/PawianParameter.hh"
#include "ErrLogger/ErrLogger.hh"
#include <vector>

class PawianParameters : public AbsPawianParameters {

public:

  PawianParameters();
  PawianParameters(std::vector<PawianParameter>& parameters);

  ~PawianParameters();

  virtual std::string type(){return "Pawian";}

  virtual unsigned int VariableParameters() const{
    return _parameters.size();
   }

  // virtual const ROOT::Minuit2::MnUserParameters& mnUserParameters() const { 
  //   MnUserParameters mnUserParameters;
  //    return  mnUserParameters;}
  
  virtual AbsPawianParameters* Clone();

  virtual AbsPawianParameters& operator=(AbsPawianParameters& par){
    if(par.type() != "Pawian"){
      Alert << "cannot copy AbsPawianParameters is type " << par.type() << " instead of type " << type() << endmsg;
      exit(1);
    }
     PawianParameters& mnPar = dynamic_cast<PawianParameters&>(par);
    _parameters = mnPar._parameters;
    return *this;
  }

  virtual PawianParameters& operator=(PawianParameters& par){
    _parameters = par._parameters;
    return *this;
  }

  virtual const std::vector<std::string> ParamNames();

   /// access to parameters and errors in column-wise representation 
  virtual std::vector<double> Params() const;
  virtual std::vector<double> Errors() const;

   /// access to single Parameter
   // const MinuitParameter& Parameter(unsigned int) const;

   /// Add free Parameter Name, Value, Error
  virtual bool Add(const std::string &, double, double);
   /// Add limited Parameter Name, Value, Lower bound, Upper bound
   virtual bool Add(const std::string &, double, double, double, double);
   /// Add const Parameter Name, vale
   virtual bool Add(const std::string &, double);

   /// interaction via external number of Parameter
   virtual void Fix(unsigned int);
   virtual void Release(unsigned int);
   virtual void SetValue(unsigned int, double);
   virtual void SetError(unsigned int, double);
   virtual void SetLimits(unsigned int, double, double);
   virtual void SetUpperLimit(unsigned int, double);
   virtual void SetLowerLimit(unsigned int, double);
   virtual void RemoveLimits(unsigned int);

   virtual double Value(unsigned int) const;
   virtual double Error(unsigned int) const;
  
   /// interaction via Name of Parameter
   virtual void Fix(const std::string &);
   virtual void Release(const std::string &);
   virtual void SetValue(const std::string &, double);
   virtual void SetError(const std::string &, double);
   virtual void SetLimits(const std::string &, double, double);
   virtual void SetUpperLimit(const std::string &, double);
   virtual void SetLowerLimit(const std::string &, double);
   virtual void RemoveLimits(const std::string &);

   virtual double Value(const std::string &) const;
   virtual double Error(const std::string &) const;

  // interaction via external number of Parameter
  virtual bool IsFixed(unsigned int);
  virtual bool HasLimits(unsigned int);
  virtual double UpperLimit(unsigned int);
  virtual double LowerLimit(unsigned int);  
  
  //convert Name into external number of Parameter
   virtual unsigned int Index(const std::string &) const;
   //convert external number into Name of Parameter
   virtual const std::string & GetName(unsigned int) const;
   // mantain interface with const char * for backward compatibility
   virtual const char* Name(unsigned int) const;

private:

  std::vector<PawianParameter> _parameters;
};

