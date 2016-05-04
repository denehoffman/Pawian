//************************************************************************//
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

#include <complex>
#include <vector>
#include "PwaDynamics/TMatrixBase.hh"

class TH2F;

class RiemannSheetAnalyzer
{

public:
   RiemannSheetAnalyzer(unsigned int noOfChannels, std::shared_ptr<TMatrixBase> tMatrix, 
			std::complex<double> massMin, std::complex<double> massMax, int numSteps, unsigned int projectionIndex=0);

private:
   void FindPolePositions(TH2F* histo);
   std::vector<double> GetSignsFromInteger(unsigned int noOfChannels, unsigned int signCollection);
};
