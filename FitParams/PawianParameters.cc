//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

#include <math.h>
#include <stdio.h>
#include <algorithm>

#include "FitParams/PawianParameters.hh"

PawianParameters::PawianParameters() :
  AbsPawianParameters()
{
}

PawianParameters::PawianParameters(std::vector<PawianParameter>& parameters) :
  AbsPawianParameters()
{
  std::vector<PawianParameter>::const_iterator itPar;
  for (itPar=parameters.begin(); itPar!=parameters.end(); ++itPar){
    PawianParameter newParam= (*itPar);
    _parameters.push_back(newParam);
  }
}

PawianParameters::~PawianParameters()
{
}

unsigned int PawianParameters::VariableParameters() const{
  unsigned int result=0;
  std::vector<PawianParameter>::const_iterator itPar;
  for (itPar=_parameters.begin(); itPar!=_parameters.end(); ++itPar){
    if (!(itPar->IsFixed())) ++result;
  }
  return result;
}

AbsPawianParameters* PawianParameters::Clone(){
  PawianParameters* result = new PawianParameters(_parameters);
  std::map<unsigned int, std::vector<unsigned int> >::iterator it;
  for(it=_depMap.begin(); it!=_depMap.end(); ++it){
    result->SetDependencies(it->first, it->second);
  }
  return result;
}

const std::vector<std::string> PawianParameters::ParamNames(){
  std::vector<std::string> parNames;
  std::vector<PawianParameter>::const_iterator itPar;
  for (itPar=_parameters.begin(); itPar!=_parameters.end(); ++itPar){
    parNames.push_back(itPar->GetName());
  }

  return parNames; 
}

bool PawianParameters::Add(const std::string & name, double val, double err){
  bool alreadyAdded=false;
  std::vector<PawianParameter>::const_iterator itPar;
  for (itPar=_parameters.begin(); itPar!=_parameters.end(); ++itPar){
    if(itPar->Name()==name) alreadyAdded=true;
  }
  if(alreadyAdded) return false; 

  _parameters.push_back(PawianParameter(_parameters.size(), name, val, err));
  return true;
}

bool PawianParameters::Add(const std::string & name, double val, double err, double low, double up){
  bool alreadyAdded=false;
  std::vector<PawianParameter>::const_iterator itPar;
  for (itPar=_parameters.begin(); itPar!=_parameters.end(); ++itPar){
    if(itPar->Name()==name) alreadyAdded=true;
  }
  if(alreadyAdded) return false; 

  _parameters.push_back(PawianParameter(_parameters.size(), name, val, err, low, up));
  return true;
}

bool PawianParameters::Add(const std::string & name, double val){
  bool alreadyAdded=false;
  std::vector<PawianParameter>::const_iterator itPar;
  for (itPar=_parameters.begin(); itPar!=_parameters.end(); ++itPar){
    if(itPar->Name()==name) alreadyAdded=true;
  }
  if(alreadyAdded) return false; 

  _parameters.push_back(PawianParameter(_parameters.size(), name, val));
  return true;
}

void PawianParameters::Fix(unsigned int idx){
  assert(idx < _parameters.size());
  _parameters[idx].Fix();
}

void PawianParameters::Release(unsigned int idx){
  assert(idx < _parameters.size());
  _parameters[idx].Release();
}

void PawianParameters::SetValue(unsigned int idx, double val){
  assert(idx < _parameters.size());
  _parameters[idx].SetValue(val);
}

void PawianParameters::SetError(unsigned int idx, double err){
  assert(idx < _parameters.size());
  _parameters[idx].SetError(err);
}

void PawianParameters::SetLimits(unsigned int idx, double low, double up){
  assert(idx < _parameters.size());
  _parameters[idx].SetLimits(low, up);
}

void PawianParameters::SetUpperLimit(unsigned int idx, double up){
  assert(idx < _parameters.size());
  _parameters[idx].SetUpperLimit(up);
}

void PawianParameters::SetLowerLimit(unsigned int idx, double low){
  assert(idx < _parameters.size());
  _parameters[idx].SetLowerLimit(low);
}

void PawianParameters::RemoveLimits(unsigned int idx){
  assert(idx < _parameters.size());
  _parameters[idx].RemoveLimits();
}

double PawianParameters::Value(unsigned int idx) const{
  assert(idx < _parameters.size());
  return _parameters[idx].Value();
}

double PawianParameters::Error(unsigned int idx) const{
  assert(idx < _parameters.size());
  return _parameters[idx].Error();
}

void PawianParameters::Fix(const std::string & pname){
  Fix(Index(pname));
}

void PawianParameters::Release(const std::string & pname){
  Release(Index(pname));
}

void PawianParameters::SetValue(const std::string & pname, double val){
  SetValue(Index(pname), val);
}

void  PawianParameters::SetError(const std::string & pname, double err){
  SetError(Index(pname), err);
}

void  PawianParameters::SetLimits(const std::string & pname, double low, double up){
  SetLimits(Index(pname), low, up);
}

void  PawianParameters::SetUpperLimit(const std::string & pname, double up){
  SetUpperLimit(Index(pname), up);
}

void  PawianParameters::SetLowerLimit(const std::string & pname, double low){
  SetLowerLimit(Index(pname), low);
}

void  PawianParameters::RemoveLimits(const std::string & pname){
  RemoveLimits(Index(pname));
}

double PawianParameters::Value(const std::string & pname) const{
  return Value(Index(pname));
}

double PawianParameters::Error(const std::string & pname) const{
  return Error(Index(pname));
}

bool PawianParameters::IsFixed(unsigned int idx){
  assert(idx < _parameters.size());
  return _parameters[idx].IsFixed();
}  

bool PawianParameters::HasLimits(unsigned int idx){
  assert(idx < _parameters.size());
  return _parameters[idx].HasLimits();
} 

double PawianParameters::UpperLimit(unsigned int idx){
  assert(idx < _parameters.size());
  return _parameters[idx].UpperLimit();
}  

double PawianParameters::LowerLimit(unsigned int idx){
  assert(idx < _parameters.size());
  return _parameters[idx].LowerLimit();
}
unsigned int PawianParameters::Index(const std::string & name) const{

  std::vector<PawianParameter>::const_iterator itPar;
  for (itPar=_parameters.begin(); itPar!=_parameters.end(); ++itPar){
    if(itPar->Name()==name) break;
  }
  assert(itPar != _parameters.end());
  return (*itPar).Number();
}

const std::string & PawianParameters::GetName(unsigned int idx) const{
  assert(idx < _parameters.size());
  return _parameters[idx].GetName();
}

const char* PawianParameters::Name(unsigned int idx) const{
  assert(idx < _parameters.size());
  return _parameters[idx].Name();
}

std::vector<double> PawianParameters::Params() const{
  // return std::vector of double with parameter values 
   unsigned int n = _parameters.size(); 
   std::vector<double> result(n);
   for(unsigned int i = 0; i < n; ++i) 
      result[i] = _parameters[i].Value(); 
   
   return result;
}

std::vector<double> PawianParameters::Errors() const{
   // return std::vector of double with parameter errors
  unsigned int n = _parameters.size(); 
   std::vector<double> result(n);
   for(unsigned int i = 0; i < n; ++i) 
      result[i] = _parameters[i].Error(); 

   return result;
}

