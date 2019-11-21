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
#include <complex>
#include "KMatrixExtract/ArgandModuloConverter.hh"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TGraphErrors.h"

#include "ErrLogger/ErrLogger.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"


ArgandModuloConverter::ArgandModuloConverter(std::string rootFileInName, 
					     double mass1, double mass2, 
					     double mass3, double mass4) :
  _theTFileOut(0)
  ,_theTFileIn(0)
  ,_php1(new PhaseSpaceIsobar(mass1, mass2))
  ,_php2(new PhaseSpaceIsobar(mass3, mass4))
{
  _theTFileIn=new TFile(rootFileInName.c_str());
  _dataTree = (TTree*) _theTFileIn->Get("_dataFourvecs");
  _fitTree = (TTree*) _theTFileIn->Get("_fittedFourvecs");

   std::string rootFileOutName="./ArgandModuloConverter.root";
   _theTFileOut=new TFile(rootFileOutName.c_str(),"recreate");

  _argandUnitsDataGraphErr = new TGraphErrors();
  _argandUnitsDataGraphErr->SetName("ArgandUnitsData");
  _argandUnitsDataGraphErr->SetTitle("ArgandUnits (data)");

  _moduloDataGraphErr = new TGraphErrors();
  _moduloDataGraphErr->SetName("ModuloData");
  _moduloDataGraphErr->SetTitle("Modulo (data)");

  _argandUnitsFitGraphErr = new TGraphErrors();
  _argandUnitsFitGraphErr->SetName("ArgandUnitsFit");
  _argandUnitsFitGraphErr->SetTitle("ArgandUnits (fit)");

  _moduloFitGraphErr = new TGraphErrors();
  _moduloFitGraphErr->SetName("ModuloFit");
  _moduloFitGraphErr->SetTitle("Modulo (fit)");
}

ArgandModuloConverter::~ArgandModuloConverter()
{
  _argandUnitsDataGraphErr->Write();
  _moduloDataGraphErr->Write();
  _argandUnitsFitGraphErr->Write();
  _moduloFitGraphErr->Write();

   _theTFileOut->Write();
   _theTFileOut->Close();
}

void ArgandModuloConverter::convertArgandToModulo(){
  _dataTree->SetBranchAddress("mass", &_massVal);
  _dataTree->SetBranchAddress("dataVal", &_dataVal);
  _dataTree->SetBranchAddress("dataErrVal", &_dataErrVal);

  int noOfEntries=_dataTree->GetEntries();
  for (int i=0; i<noOfEntries; ++i){
    _dataTree->GetEntry(i);
    _argandUnitsDataGraphErr->SetPoint(i, _massVal, _dataVal);
    _argandUnitsDataGraphErr->SetPointError(i, 0., _dataErrVal);

    double factor=1./norm(sqrt(_php1->factor(_massVal)*_php2->factor(_massVal)));
    _moduloDataVal=factor*_dataVal;
    _moduloDataErrVal=factor*_dataErrVal;
    _moduloDataGraphErr->SetPoint(i, _massVal, _moduloDataVal);
    _moduloDataGraphErr->SetPointError(i, 0., _moduloDataErrVal);
  }

  _fitTree->SetBranchAddress("mass", &_massVal);
  _fitTree->SetBranchAddress("fitVal", &_fitVal);

  noOfEntries=_fitTree->GetEntries();
  for (int i=0; i<noOfEntries; ++i){
    _fitTree->GetEntry(i);
    _argandUnitsFitGraphErr->SetPoint(i, _massVal, _fitVal);
    _argandUnitsFitGraphErr->SetPointError(i, 0., 0.);
    double factor=1./norm(sqrt(_php1->factor(_massVal)*_php2->factor(_massVal)));
    _moduloFitVal=factor*_fitVal;
    _moduloFitGraphErr->SetPoint(i, _massVal, _moduloFitVal);
    _moduloFitGraphErr->SetPointError(i, 0., 0.);
  }
}

