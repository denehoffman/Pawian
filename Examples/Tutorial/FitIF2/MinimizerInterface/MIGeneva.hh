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

#ifndef _MIGENEVA_H
#define _MIGENEVA_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include "Examples/Tutorial/FitIF2/MinimizerInterface/MIGeneric.hh"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// Geneva header files go here
#include <hap/GRandomT.hpp>
#include <common/GCommonEnums.hpp>
#include <common/GExceptions.hpp>
#include <geneva/GConstrainedDoubleObject.hpp>
#include <geneva/GConstrainedDoubleObjectCollection.hpp>
#include <geneva/GDoubleGaussAdaptor.hpp>
#include <geneva/GObjectExpectationChecksT.hpp>
#include <geneva/GParameterObjectCollection.hpp>
#include <geneva/GParameterSet.hpp>
#include <courtier/GAsioHelperFunctions.hpp>
#include <courtier/GAsioTCPClientT.hpp>
#include <courtier/GAsioTCPConsumerT.hpp>
#include <geneva/GParameterSet.hpp>
#include <geneva/GBrokerEA.hpp>
#include <geneva/GEvolutionaryAlgorithm.hpp>
#include <geneva/GIndividual.hpp>
#include <geneva/GMultiThreadedEA.hpp>
#include <common/GSerializationHelperFunctionsT.hpp>
#include <geneva/GOptimizationEnums.hpp>

using namespace std;

class MIGeneva : public Gem::Geneva::GParameterSet {

  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) {
    using boost::serialization::make_nvp;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GParameterSet);

		/* Add your own class-variables here in the following way:
			ar & BOOST_SERIALIZATION_NVP(myVar);
			or
			ar & make_nvp("myVar", myVar); // The latter form can be necessary when dealing with templates
		 */
  }

public:
  // Constructs the object from a generic evaluator
  MIGeneva(boost::shared_ptr<MIGeneric> theData);  

  // Gives access to a GenericInterface object from the converged solution
  const MIGeneric& getResult() const;

  /** Destructor */
  virtual ~MIGeneva();

 protected:

  virtual double fitnessCalculation();

 private:
  boost::shared_ptr<MIGeneric> genericEvaluator_;  
 // boost::shared_ptr<MIData> _myData;
 // vector<string> paramNames;

};

#endif /* _MIGeneva_H */
