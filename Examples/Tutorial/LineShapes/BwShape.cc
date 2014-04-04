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
#include "Examples/Tutorial/LineShapes/BwShape.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"
#include "PwaDynamics/BreitWignerFunction.hh"

BwShape::BwShape(double MassRes, double MassWidth, double MassDec1, double MassDec2, unsigned int Lmax, double deltaMass) :
  _theTFile(0)
{

  if (MassRes <= MassDec1+MassDec2){
    Alert << "Mass of the resonance smaller than the masses of the decay particles\n" 
	  << MassRes << " < " << MassDec1+MassDec2 << endmsg;
    exit(0); 
  }

   std::stringstream Lmaxstrstr;
   Lmaxstrstr << Lmax;
   std::string rootFileName="./BwShapeLmax"+Lmaxstrstr.str()+".root";
   _theTFile=new TFile(rootFileName.c_str(),"recreate");
 
  for (unsigned int i=0; i<=Lmax; ++i){

    std::stringstream Lstrstr;
    Lstrstr << i;
    std::string histName="BreitWigner_L"+Lstrstr.str();


   _histMap[i]= new TH1F(histName.c_str(),histName.c_str(),301, MassRes-deltaMass, MassRes+deltaMass);

   std::string histNameNew="BreitWigner_Lnew"+Lstrstr.str();
   _histMapNew[i]= new TH1F(histNameNew.c_str(),histNameNew.c_str(),301, MassRes-deltaMass, MassRes+deltaMass);

    histName="Argand_L"+Lstrstr.str();
    TH2F* currentArgandHist=new TH2F(histName.c_str(),histName.c_str(),301, -1., 1., 301, 0., 1.3);
    currentArgandHist->SetXTitle("Re(Bw)");
    currentArgandHist->SetYTitle("Im(Bw)");
    currentArgandHist->SetMarkerStyle(6);
    _argandHistMap[i]=currentArgandHist;
    
    histNameNew="Argand_Lnew"+Lstrstr.str();
    TH2F* currentArgandHistNew=new TH2F(histNameNew.c_str(),histNameNew.c_str(),301, -1., 1., 301, 0., 1.3);
    currentArgandHistNew->SetXTitle("Re(Bw)");
    currentArgandHistNew->SetYTitle("Im(Bw)");
    currentArgandHistNew->SetMarkerStyle(6);
    _argandHistMapNew[i]=currentArgandHistNew;
    
    //initialize here the phase shift histogramms _phaseHistMap
  }
  
  double stepSize=2.*deltaMass/300;
  
  for (unsigned int lIt=0; lIt<=Lmax; ++lIt){
    
    TH1F* currentHist=_histMap[lIt];
    TH1F* currentHistNew=_histMapNew[lIt];
    TH2F* currentArgandHist=_argandHistMap[lIt];

    TH2F* currentArgandHistNew=_argandHistMapNew[lIt];

    for (double massIt=MassRes-deltaMass; massIt<MassRes+deltaMass; massIt+=stepSize){
      
      Vector4<double>  res4V(massIt, 0., 0., 0.);
      
      complex<double> currentBW=BreitWignerBlattW(res4V, MassDec1, MassDec2, MassRes, MassWidth, lIt);
      complex<double> currentBWnew=BreitWignerFunction::BlattWRel(lIt, res4V.Mass(), MassRes, MassWidth, MassDec1, MassDec2);
      double weight=norm(currentBW);
      currentHist->Fill(massIt,weight);
      currentHistNew->Fill(massIt,norm(currentBWnew));
      
      currentArgandHist->Fill(currentBW.real(),currentBW.imag());

      currentArgandHistNew->Fill(currentBWnew.real(),currentBWnew.imag());
      
      //fill here the phase shift histogramms: hint "atan2(imag,real)"

    }

  }
}

BwShape::~BwShape()
{
   _theTFile->Write();
   _theTFile->Close();
}

