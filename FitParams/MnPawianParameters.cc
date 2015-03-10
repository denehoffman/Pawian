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
#include <boost/timer/timer.hpp>

#include "FitParams/MnPawianParameters.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MinuitParameter.h"

using namespace ROOT::Minuit2;

MnPawianParameters::MnPawianParameters() :
  AbsPawianParameters()
  ,_mnUserParameters()
{
}

MnPawianParameters::MnPawianParameters(const MnUserParameters& mnUserParameters) :
  AbsPawianParameters()
  ,_mnUserParameters(mnUserParameters)
{
}

MnPawianParameters::~MnPawianParameters()
{
}


AbsPawianParameters* MnPawianParameters::Clone(){
  MnPawianParameters* result = new MnPawianParameters(_mnUserParameters); 
  return result;
}


const std::vector<std::string> MnPawianParameters::ParamNames(){
  const std::vector<MinuitParameter> theParams= _mnUserParameters.Parameters();
  std::vector<std::string> parNames;
  std::vector<MinuitParameter>::const_iterator itPar;
  for (itPar=theParams.begin(); itPar!=theParams.end(); ++itPar){
    parNames.push_back(itPar->GetName());
  }

  return parNames; 
}

bool MnPawianParameters::Add(const std::string & name, double val, double err){
  return _mnUserParameters.Add(name, val, err);
}

bool MnPawianParameters::Add(const std::string & name, double val, double err, double low, double up){
  return _mnUserParameters.Add(name, val, err, low, up);
}

bool MnPawianParameters::Add(const std::string & name, double val){
  return _mnUserParameters.Add(name, val);
}

void MnPawianParameters::Fix(unsigned int idx){
  _mnUserParameters.Fix(idx);
}

void MnPawianParameters::Release(unsigned int idx){
  _mnUserParameters.Release(idx);
}

void MnPawianParameters::SetValue(unsigned int idx, double val){
  _mnUserParameters.SetValue(idx, val);
}

void MnPawianParameters::SetError(unsigned int idx, double err){
  _mnUserParameters.SetValue(idx, err);
}

void MnPawianParameters::SetLimits(unsigned int idx, double low, double up){
  _mnUserParameters.SetLimits(idx, low, up);
}

void MnPawianParameters::SetUpperLimit(unsigned int idx, double up){
  _mnUserParameters.SetUpperLimit(idx, up);
}

void MnPawianParameters::SetLowerLimit(unsigned int idx, double low){
  _mnUserParameters.SetLowerLimit(idx, low);
}

void MnPawianParameters::RemoveLimits(unsigned int idx){
  _mnUserParameters.RemoveLimits(idx);
}

double MnPawianParameters::Value(unsigned int idx) const{
  return _mnUserParameters.Value(idx);
}

double MnPawianParameters::Error(unsigned int idx) const{
  return _mnUserParameters.Error(idx);
}

void MnPawianParameters::Fix(const std::string & pname){
  _mnUserParameters.Fix(pname);
}

void MnPawianParameters::Release(const std::string & pname){
  _mnUserParameters.Release(pname);
}

void MnPawianParameters::SetValue(const std::string & pname, double val){
  _mnUserParameters.SetValue(pname, val);
}

void  MnPawianParameters::SetError(const std::string & pname, double err){
 _mnUserParameters.SetValue(pname, err);
}

void  MnPawianParameters::SetLimits(const std::string & pname, double low, double up){
  _mnUserParameters.SetLimits(pname, low, up);
}

void  MnPawianParameters::SetUpperLimit(const std::string & pname, double up){
  _mnUserParameters.SetUpperLimit(pname, up);
}

void  MnPawianParameters::SetLowerLimit(const std::string & pname, double low){
  _mnUserParameters.SetLowerLimit(pname, low);
}

void  MnPawianParameters::RemoveLimits(const std::string & pname){
  _mnUserParameters.RemoveLimits(pname);
}

double MnPawianParameters::Value(const std::string & pname) const{
  return _mnUserParameters.Value(pname);
}

double MnPawianParameters::Error(const std::string & pname) const{
  return _mnUserParameters.Error(pname);
}

bool MnPawianParameters::IsFixed(unsigned int idx){
  return _mnUserParameters.Parameter(idx).IsFixed();
}  

bool MnPawianParameters::HasLimits(unsigned int idx){
  return _mnUserParameters.Parameter(idx).HasLimits();
} 

double MnPawianParameters::UpperLimit(unsigned int idx){
  return _mnUserParameters.Parameter(idx).UpperLimit();
}  

double MnPawianParameters::LowerLimit(unsigned int idx){
  return _mnUserParameters.Parameter(idx).LowerLimit();
}
unsigned int MnPawianParameters::Index(const std::string & name) const{
  return _mnUserParameters.Index(name);
}

const std::string & MnPawianParameters::GetName(unsigned int idx) const{
  return _mnUserParameters.GetName(idx);
}

const char* MnPawianParameters::Name(unsigned int idx) const{
  return _mnUserParameters.Name(idx);
}

void MnPawianParameters::SetPrecision(double eps){
  _mnUserParameters.SetPrecision(eps);
}
