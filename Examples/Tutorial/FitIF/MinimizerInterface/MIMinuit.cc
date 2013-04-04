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

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIMinuit.hh"
//#include "Examples/Tutorial/FitIF/MinimizerInterface/MIBase.hh"

using namespace ROOT::Minuit2;
using namespace std;

MIMinuit::MIMinuit(boost::shared_ptr<MIData> theData) : _myFcn(theData){
  //_myFcn = new MIMinuitFcn(theData);
}

MIMinuit::~MIMinuit()
{
  //delete _myFcn;
}

const double MIMinuit::exec(unsigned int Npar,  double* par,  double* min, double* max, double* err){

  string s;
  stringstream out;
  
  MnUserParameters upar;
  for(unsigned int i=0; i<Npar; ++i){
    out.str("");
    out << i;
    s = out.str();
    upar.Add(s, par[i], err[i], max[i], min[i]);
  }

  MnMigrad migrad(_myFcn, upar);
 // Info <<"start migrad "<< endmsg;
  FunctionMinimum minMin = migrad();

 if(!minMin.IsValid()) {
   //try with higher strategy
 //  Info <<"FM is invalid, try with strategy = 2."<< endmsg;
   MnMigrad migrad2(_myFcn, minMin.UserState(), MnStrategy(2));
   minMin = migrad2();
 }
 
  //save minized values
  for(unsigned int i=0; i<Npar; ++i){
    out.str("");
    out << i;
    s = out.str();
    par[i]=minMin.UserState().Value(s); err[i]=minMin.UserState().Error(s);
  }

  return minMin.Fval();
}
