//**************************************************************************/
//									  //
//  Copyright 2021 Bertram Kopf (bertram@ep1.rub.de)			  //
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


#include "KMatrixExtract/RiemannSheetFVectorAnalyzer.hh"
#include "TH2F.h"
#include "ErrLogger/ErrLogger.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"


RiemannSheetFVectorAnalyzer::RiemannSheetFVectorAnalyzer(unsigned int noOfChannels, 
							 std::shared_ptr<FVector> fVector,
							 std::complex<double> massMin, 
							 std::complex<double> massMax,
							 int numXSteps,
							 int numYSteps,
							 unsigned int projectionIndex)
{
   unsigned int signCollection = 0;

   for(unsigned int i=0; i<pow(2, noOfChannels); i++) { 
      std::vector<double> signs = GetSignsFromInteger(noOfChannels, signCollection);
      std::ostringstream histoname;
    
      histoname << "sheet_";
      // bool fillSheet=true;
      // bool lastSignIsNeg=true;
     
      for(auto it = signs.begin(); it!= signs.end();++it){
	 histoname << ((*it) == 1 ? "+" : "-");
	 // bool currentSignIsNeg=true;
         // if( (*it) > 0. ) currentSignIsNeg=false;
         // if ( (!lastSignIsNeg) && currentSignIsNeg ) fillSheet=false;
	 
         // lastSignIsNeg=currentSignIsNeg;
	 
      // }

      // if( !fillSheet ){
      // 	signCollection++;
      // 	continue;
      }

      InfoMsg << "Scanning " << histoname.str() << endmsg;      
      TH2F* scan = new TH2F(histoname.str().c_str(), histoname.str().c_str(), 
			    numXSteps, massMin.real(), massMax.real(), 
			    numYSteps, massMin.imag(), massMax.imag());

      fVector->SetBumImPartSigns(signs);

      for(int i=1;i<=scan->GetNbinsX(); i++){
	 for(int j=1;j<=scan->GetNbinsY(); j++){
	    fVector->evalMatrix(std::complex<double>(scan->GetXaxis()->GetBinCenter(i), 
						     scan->GetYaxis()->GetBinCenter(j)));
	    // scan->SetBinContent(i,j, std::abs((*tMatrix)(projectionIndex,projectionIndex)));
	    complex<double> currentrho
	      = fVector->kMatrix()->phaseSpaceVec().at(projectionIndex)->factor(std::complex<double>(scan->GetXaxis()->GetBinCenter(i), scan->GetYaxis()->GetBinCenter(j)));
	    scan->SetBinContent(i,j, std::abs((*fVector)(projectionIndex,0)
					      *sqrt(currentrho)));
	 }
      }

      FindPolePositions(scan);

      InfoMsg << "\n" << endmsg;

      signCollection++;
   }
}


