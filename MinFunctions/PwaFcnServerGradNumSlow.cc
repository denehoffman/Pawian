//************************************************************************//
//									  //
//  Copyright 2024 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include <iomanip>
#include "Minuit2/FCNBase.h"
#include "Minuit2/FCNGradientBase.h"

#include "MinFunctions/PwaFcnServerGradNumSlow.hh"
#include "FitParams/ParamDepHandler.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/NetworkServer.hh"
#include "Utils/PawianConstants.hh"
#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

PwaFcnServerGradNumSlow::PwaFcnServerGradNumSlow(std::shared_ptr<NetworkServer> netServer) :
  PwaFcnServer<FCNGradientBase>(netServer)
{
}

PwaFcnServerGradNumSlow::~PwaFcnServerGradNumSlow()
{
}

std::vector<double> PwaFcnServerGradNumSlow::Gradient(const std::vector<double>& par) const{
  std::vector<double> resultVec;
  resultVec.resize(par.size());
  this->_currentPawianParms->SetAllValues(par);
  ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);
  double LHBase=collectLH();
  double epsilon;
  
  for(unsigned int i=0; i<par.size(); ++i){
    if(this->_currentPawianParms->IsFixed(i)) resultVec.at(i)=0.;
    else{
      double currentVal=this->_currentPawianParms->Value(i);
      epsilon=_numStepSize*std::abs(currentVal);
      if ((this->_currentPawianParms->GetName(i)).substr( (this->_currentPawianParms->GetName(i)).length() - 3 ) == "Phi"){
	epsilon=_numStepSize*PawianConstants::pi;
      }
      else if (std::abs(currentVal)<1.e-10)  epsilon=_numStepSize*1.e-10;
      double dx=(currentVal+epsilon)-currentVal;
      
      double Lh_high=0.;
      double Lh_low=0.;
      if(!(this->_currentPawianParms->HasLimits(i)) || (std::abs(currentVal-this->_currentPawianParms->UpperLimit(i))>epsilon  && std::abs(currentVal-this->_currentPawianParms->LowerLimit(i))>epsilon) ){
	this->_currentPawianParms->SetValue(i, currentVal+epsilon);
	ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);
	Lh_high=collectLH();
	this->_currentPawianParms->SetValue(i, currentVal-epsilon);
	ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);
	Lh_low=collectLH();
	resultVec.at(i)=(Lh_high-Lh_low)/(2.*dx);
      }
      else{
	if(std::abs(currentVal-this->_currentPawianParms->UpperLimit(i))<epsilon){
	  this->_currentPawianParms->SetValue(i, currentVal-epsilon);
	  ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);
	  Lh_low=collectLH();
	  resultVec.at(i)=(LHBase-Lh_low)/dx;
	}
	else{
	  this->_currentPawianParms->SetValue(i, currentVal+epsilon);
	  ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);
	  Lh_high=collectLH();
	  resultVec.at(i)=(Lh_high-LHBase)/dx;
	}
      }
      this->_currentPawianParms->SetValue(i, currentVal);      
    }
  }
  return resultVec;  
}





