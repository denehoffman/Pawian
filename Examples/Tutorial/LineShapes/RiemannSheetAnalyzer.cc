/**************************************************************************/
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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


#include "Examples/Tutorial/LineShapes/RiemannSheetAnalyzer.hh"
#include "TH2F.h"
#include "ErrLogger/ErrLogger.hh"

const unsigned int NUM_BINS = 500;

RiemannSheetAnalyzer::RiemannSheetAnalyzer(unsigned int noOfChannels, 
					   std::shared_ptr<TMatrixBase> tMatrix,
					   std::complex<double> massMin, 
					   std::complex<double> massMax)
{
   unsigned int signCollection = 0;

   for(unsigned int i=0; i<pow(2, noOfChannels); i++){ 
      std::vector<double> signs = GetSignsFromInteger(noOfChannels, signCollection);
      std::ostringstream histoname;
    
      histoname << "sheet_";
      for(auto it = signs.begin(); it!= signs.end();++it){
	 histoname << ((*it) == 1 ? "+" : "-");
      }
      
      Info << "Scanning " << histoname.str() << endmsg;

      TH2F* scan = new TH2F(histoname.str().c_str(), histoname.str().c_str(), 
			    NUM_BINS, massMin.real(), massMax.real(), 
			    NUM_BINS, massMin.imag(), massMax.imag());

      tMatrix->SetBumImPartSigns(signs);
      
      for(int i=1;i<=scan->GetNbinsX(); i++){
	 for(int j=1;j<=scan->GetNbinsY(); j++){
	    tMatrix->evalMatrix(std::complex<double>(scan->GetXaxis()->GetBinCenter(i), scan->GetYaxis()->GetBinCenter(j)));
	    scan->SetBinContent(i,j, std::abs((*tMatrix)(0,0)));		 
	 }
      }

      FindPolePositions(scan);

      Info << "\n" << endmsg;

      signCollection++;
   }
}



std::vector<double> RiemannSheetAnalyzer::GetSignsFromInteger(unsigned int noOfChannels, unsigned int signCollection){

   std::vector<double> result;

   for(unsigned int i=0; i<noOfChannels;i++){
      unsigned int currentBit = (signCollection >> i) & 0x1;
      if(currentBit == 1){
	 result.push_back(1);
      }
      else{
	 result.push_back(-1);
      }
   }

   return result;
}



void RiemannSheetAnalyzer::FindPolePositions(TH2F* histo){

   unsigned int numPolesFound=0;

   for(int i=1;i<=histo->GetNbinsX(); i++){
      for(int j=1;j<=histo->GetNbinsY(); j++){
	 if(i==1 || j==1 || i==histo->GetNbinsX() || j==histo->GetNbinsY())
	    continue;

	 double val = histo->GetBinContent(i, j);
	 if( (val > histo->GetBinContent(i+1, j)) && (val > histo->GetBinContent(i-1, j)) &&
	     (val > histo->GetBinContent(i, j+1)) && (val > histo->GetBinContent(i, j-1))){
	    Info << "Found possible pole at (" 
		 << histo->GetXaxis()->GetBinCenter(i) << ", "
		 << -1*2*histo->GetYaxis()->GetBinCenter(j) << ") GeV" << endmsg;
	    numPolesFound++;
	 }
      }
   }

   if(numPolesFound > 0){
      Info << "[binning uncertainty (" << histo->GetXaxis()->GetBinWidth(1) << ", "
	   << histo->GetYaxis()->GetBinUpEdge(1) - histo->GetYaxis()->GetBinLowEdge(1) << ")]" << endmsg;
   }
}
