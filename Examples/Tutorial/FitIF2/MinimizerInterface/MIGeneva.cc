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

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "Examples/Tutorial/FitIF2/MinimizerInterface/MIGeneva.hh"
#include "Examples/Tutorial/FitIF2/MinimizerInterface/MIGeneric.hh"

using namespace std;
using namespace Gem::Geneva;
using namespace Gem::Courtier;
using namespace Gem::Hap;

MIGeneva::MIGeneva(boost::shared_ptr<MIGeneric> theData): GParameterSet(){

  genericEvaluator_=theData;

  // Add bounded double objects
  for(std::size_t i=0; i<genericEvaluator_->getNumPars(); i++) {
			
    boost::shared_ptr<GConstrainedDoubleObject> gbd_ptr(new GConstrainedDoubleObject(genericEvaluator_->getPar(i), genericEvaluator_->getMin(i), genericEvaluator_->getMax(i)) );

    // Create a suitable adaptor (sigma=0.1, sigma-adaption=0.5, min sigma=0, max sigma=0,5)
    boost::shared_ptr<GDoubleGaussAdaptor> gdga_ptr(new GDoubleGaussAdaptor(genericEvaluator_->getErr(i), 0.5, 0., 3*genericEvaluator_->getErr(i)));
    gdga_ptr->setAdaptionThreshold(1); // Adaption parameters are modified after each adaption
    gdga_ptr->setAdaptionProbability(0.05); // The likelihood for a parameter to be adapted

    // Register the adaptor with GConstrainedDoubleObject objects
    gbd_ptr->addAdaptor(gdga_ptr);

    // Add a GConstrainedDoubleObject object to the collection
    // gbdc_ptr->push_back(gbd_ptr);
    // gpoc_ptr->push_back(gbd_ptr);
    this->push_back(gbd_ptr);
  }

}

MIGeneva::~MIGeneva()
{
  //delete _myFcn;
}

double MIGeneva::fitnessCalculation(){
  //double result = 0.;

  // Extract the GConstrainedDoubleObjectCollection object. In a realistic scenario, you might want
  // to add error checks here upon first invocation.
      
  return genericEvaluator_->evaluator();
}

const MIGeneric& MIGeneva::getResult() const {


  return *genericEvaluator_;// bestIndividual_ptr->fitnessCalculation();
}
