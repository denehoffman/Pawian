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

#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Tutorial/LineShapes/VoigtShape.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
//#include "TMath.h"
//#include "RooVoigtian.h"
//#include "RooAbsReal.h"
//#include "RooRealVar.h"

#include "ErrLogger/ErrLogger.hh"
#include "PwaDynamics/Voigtian.hh"

VoigtShape::VoigtShape(double mass0, double gamma, double sigma) :
  _theTFile(0)
{


   std::string rootFileName="./Voigtian.root";
   _theTFile=new TFile(rootFileName.c_str(),"recreate");
   double deltaMass=6.*(gamma+sigma);
   if(mass0-deltaMass<0.) deltaMass=mass0;

   double massMin=mass0-deltaMass;
   double massMax=mass0+deltaMass;
   TH1F* _voigtMassH1=new TH1F("voigtMassH1", "voigtMassH1",301, massMin, massMax);
   TH1F* _bwMassH1=new TH1F("bwMassH1", "bwMassH1",301, massMin, massMax);
   Voigtian theVoigtian;  
 
  double stepSize=2.*deltaMass/300;
  
  for (double massIt=mass0-deltaMass; massIt<mass0+deltaMass; massIt+=stepSize){
    
    Vector4<double>  res4V(massIt, 0., 0., 0.);
    double currentVoigtian=theVoigtian.calc(massIt, mass0, gamma, sigma);
    InfoMsg << "currentVoigtian(" << massIt << ")= " << currentVoigtian << endmsg;    
    double weight=currentVoigtian;
    _voigtMassH1->Fill(massIt,weight);
    
    complex<double> currentBW=BreitWigner(res4V,mass0, gamma);
    double weightBW=norm(currentBW);
    _bwMassH1->Fill(massIt,weightBW);
     
  }
}

VoigtShape::~VoigtShape()
{
   _theTFile->Write();
   _theTFile->Close();
}

