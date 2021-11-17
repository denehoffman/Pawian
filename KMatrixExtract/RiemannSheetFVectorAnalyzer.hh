//************************************************************************//
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

#include <complex>
#include <vector>

#include "KMatrixExtract/RiemannSheetAnalyzer.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/FVector.hh"
class TH2F;

class RiemannSheetFVectorAnalyzer : public RiemannSheetAnalyzer
{

public:
   RiemannSheetFVectorAnalyzer(unsigned int noOfChannels, std::shared_ptr<FVector> fVector, 
			       std::complex<double> massMin, std::complex<double> massMax, int numXSteps, int numYSteps, unsigned int projectionIndex=0);

  RiemannSheetFVectorAnalyzer(){;}
  
private:
};
