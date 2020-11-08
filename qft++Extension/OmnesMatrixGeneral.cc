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

#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "qft++Extension/OmnesMatrixGeneral.hh"
#include "qft++Extension/AbsOmnesMatrix.hh"
#include "qft++Extension/AbsComplexFunction.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

OmnesMatrixGeneral::OmnesMatrixGeneral(std::shared_ptr<AbsOmnesMatrix> absOmnesMatr) :
  _absOmnesMatr(absOmnesMatr)
  ,_omnesFunctions(boost::extents[absOmnesMatr->omnesMatrix().size()][absOmnesMatr->omnesMatrix().size()])
  ,_selfEnergyFunctions(boost::extents[absOmnesMatr->omnesMatrix().size()][absOmnesMatr->omnesMatrix().size()])
  ,_theTFile(0)
{
    std::string rootFileName="./OmnesMatrixGeneral.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");
  
  int omnesDimension=absOmnesMatr->omnesMatrix().size();
  for(int i=0; i<omnesDimension; ++i){
    std::stringstream strstri;
    strstri << i;
    std::string stri = strstri.str();
    for(int j=0; j<omnesDimension; ++j){
      std::stringstream strstrj;
      strstrj << j;
      std::string strj = strstrj.str();
      std::string omnesStr="Omnes";

      std::string currentOmnesRealName=omnesStr+"Real_"+stri+"_"+strj;
      TH1F* currentOmnesRealH1=new TH1F(currentOmnesRealName.c_str(), currentOmnesRealName.c_str(), 
                                    2000-1, 0., 2.);
      _omnesRealH1Vec[i].push_back(currentOmnesRealH1);

      std::string currentOmnesImagName=omnesStr+"Imag_"+stri+"_"+strj;
      TH1F* currentOmnesImagH1=new TH1F(currentOmnesImagName.c_str(), currentOmnesImagName.c_str(),
                                    2000-1, 0., 2.);
      _omnesImagH1Vec[i].push_back(currentOmnesImagH1);

      std::string selfEnergyStr="Self-Energy";
      std::string currentSelfERealName=selfEnergyStr+"Real_"+stri+"_"+strj;
      TH1F* currentSelfERealH1=new TH1F(currentSelfERealName.c_str(), currentSelfERealName.c_str(),
                                    2000-1, 0., 2.);
      _selfEnergyRealH1Vec[i].push_back(currentSelfERealH1);

      std::string currentSelfEImagName=selfEnergyStr+"Imag_"+stri+"_"+strj;
      TH1F* currentSelfEImagH1=new TH1F(currentSelfEImagName.c_str(), currentSelfEImagName.c_str(),
                                    2000-1, 0., 2.);
      _selfEnergyImagH1Vec[i].push_back(currentSelfEImagH1);
      
      _omnesFunctions[i][j]=absOmnesMatr->omnesMatrix()[i][j];
      _selfEnergyFunctions[i][j]=absOmnesMatr->selfEnergyMatrix()[i][j];
    }
  }
  
  InfoMsg << "omnesDimension: " << omnesDimension << endmsg;
    
  for(double currentMass=0.; currentMass<2.; currentMass += 0.001){
    for(int i=0; i<omnesDimension; ++i){
      for(int j=0; j<omnesDimension; ++j){
	complex<double> currentOmnesVal=_omnesFunctions[i][j]->eval(currentMass*currentMass);
	_omnesRealH1Vec.at(i).at(j)->Fill(currentMass, currentOmnesVal.real());
	_omnesImagH1Vec.at(i).at(j)->Fill(currentMass, currentOmnesVal.imag());
	
        complex<double> currentSelfEVal=_selfEnergyFunctions[i][j]->eval(currentMass*currentMass);
        _selfEnergyRealH1Vec.at(i).at(j)->Fill(currentMass, currentSelfEVal.real());
        _selfEnergyImagH1Vec.at(i).at(j)->Fill(currentMass, currentSelfEVal.imag());

      }
    }
  }

}

OmnesMatrixGeneral::~OmnesMatrixGeneral() {
    _theTFile->Write();
    _theTFile->Close();
}


