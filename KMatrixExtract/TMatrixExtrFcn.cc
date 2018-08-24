//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>  // std::setprecision
//#include <boost/multi_array.hpp>
#include "KMatrixExtract/TMatrixExtrFcn.hh"
#include "KMatrixExtract/TMatrixExtrFit.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

TMatrixExtrFcn::TMatrixExtrFcn(std::shared_ptr<TMatrixExtrFit> tMatFit) :
  _tMatFit(tMatFit)
{ 
}

TMatrixExtrFcn::TMatrixExtrFcn(std::shared_ptr<TMatrixErrorExtr> tMatErrorExtr) :
  _tMatErrorExtr(tMatErrorExtr)
{ 
}

TMatrixExtrFcn::~TMatrixExtrFcn()
{
}

double TMatrixExtrFcn::operator()(const std::vector<double>& par) const
{
  double result=-log(_tMatFit->calcTMatrix(par.at(0), par.at(1)));
  //double result=1./_tMatFit->calcTMatrix(par.at(0), par.at(1));  
  InfoMsg << "par0: " << par.at(0) << "\tpar1: " << par.at(1) 
	  << "\tabs(tMat): " << abs(_tMatFit->calcTMatrix(par.at(0), par.at(1)))
	  << "\tresult: " << result << endmsg;
  return result;
}

double TMatrixExtrFcn::Up() const 
{
  return 1.;
}

