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
#include <map>
#include <complex>
#include "TFile.h"
#include "TH1F.h"

#include "Examples/Tutorial/DfuncClebschG/SphericalHarmonicsPolts.hh"
#include "qft++Extension/SphericalHarmonic.hh"


#include "ErrLogger/ErrLogger.hh"


SphericalHarmonicsPolts::SphericalHarmonicsPolts(int lmax){
  _theTFile=new TFile("SphericalHarmonicsPolts.root","recreate");

  for (int theL=0; theL<=lmax; ++theL){
    for (int mL=-theL; mL<=theL; ++mL){
      std::stringstream LStrStr;
      LStrStr << theL;
      std::stringstream mLStrStr;
      mLStrStr << mL;
      std::string currentHistName="Y_"+LStrStr.str()+"_"+mLStrStr.str();
      TH1F* currentHist = new TH1F(currentHistName.c_str(), currentHistName.c_str(), 201., -1., 1.);
      _YLmLHists[theL][mL]=currentHist;
    }
  }

  for (double cosT=-1.; cosT<=1.; cosT+=2./200.){
    double theta=acos(cosT);
    for (int theL=0; theL<=lmax; ++theL){
      for (int mL=-theL; mL<=theL; ++mL){
	std::complex<double> sphereHarmAmp=SphericalHarmonic::SphericalHarmonicFunc(theL, mL, theta, 0.);
	_YLmLHists.at(theL).at(mL)->Fill(cosT, norm(sphereHarmAmp));	
      }
    }
  }
}

SphericalHarmonicsPolts::~SphericalHarmonicsPolts()
{
  _theTFile->Write();
  _theTFile->Close();
}



