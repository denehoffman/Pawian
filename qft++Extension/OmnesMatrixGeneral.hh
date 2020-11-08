//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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
//OmnesMatrixGeneral class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <memory>

#include <boost/multi_array.hpp>

#include "TROOT.h"


class TFile;
class TH1F;
class TH2F;
class AbsOmnesMatrix;
class AbsComplexFunction;

class OmnesMatrixGeneral {

public:

  // create/copy/destroy:

  ///Constructor
  OmnesMatrixGeneral(std::shared_ptr<AbsOmnesMatrix> absOmnesMatr); 


  /** Destructor */
  virtual ~OmnesMatrixGeneral();

  //  virtual void process();

  // Getters:
 
protected:
  std::shared_ptr<AbsOmnesMatrix> _absOmnesMatr;
  boost::multi_array< std::shared_ptr<AbsComplexFunction> , 2> _omnesFunctions;
  boost::multi_array< std::shared_ptr<AbsComplexFunction> , 2> _selfEnergyFunctions;
  
  TFile* _theTFile;

private:
  std::map< int, std::vector<TH1F*> > _omnesRealH1Vec;
  std::map< int, std::vector<TH1F*> > _omnesImagH1Vec;
  std::map< int, std::vector<TH1F*> > _selfEnergyRealH1Vec;
  std::map< int,std::vector<TH1F*> > _selfEnergyImagH1Vec;

  void init();
};


