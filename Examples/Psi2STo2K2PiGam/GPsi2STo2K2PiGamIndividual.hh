// Standard header files go here
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm> // for std::sort
#include <utility> // For std::pair
#include <cassert>


// Includes check for correct Boost version(s)
#include "GGlobalDefines.hpp"

// Boost header files go here
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/cast.hpp>

#ifndef GPsi2STo2K2PiGamIndividual_HPP_
#define GPsi2STo2K2PiGamIndividual_HPP_

// For Microsoft-compatible compilers
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

// GenEvA header files go here
#include "GParameterSet.hpp"
#include "GDoubleGaussAdaptor.hpp"
#include "GBoundedDoubleCollection.hpp"
#include "GBoundedDouble.hpp"
#include "GenevaExceptions.hpp"
#include "GGlobalOptionsT.hpp"
#include "GPODExpectationChecksT.hpp"
#include "GEnums.hpp"
#include "GUnitTestFrameworkT.hpp"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MinuitParameter.h"

#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"

using namespace ROOT::Minuit2;

namespace Gem
{
namespace GenEvA
{
/************************************************************************************************/
/**
 * This individual searches for the minimum of a parabola of a given dimension,
 * It is part of a complete example that lets users adapt their optimization
 * problems more easily to the Geneva conventions.
 */
class GPsi2STo2K2PiGamIndividual :public GParameterSet
{
	///////////////////////////////////////////////////////////////////////
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
	///////////////////////////////////////////////////////////////////////

public:
	/********************************************************************************************/
	/**
	 * A simple constructor that initializes this object with a collection of bounded
	 * double variables.
	 *
	 * @param dim The amount of variables
	 * @param min The lower boundary of the variables
	 * @param max The upper boundary of the variables
	 */
	GPsi2STo2K2PiGamIndividual(boost::shared_ptr<const AbsPsi2STo2K2PiGamLh> theLh, MnUserParameters& upar)
	  : GParameterSet()
	  ,_psi2STo2K2PiGamLhPtr(theLh->clone_())
	  ,_upar(upar)
	  {
	    setUp();
	  }
  
  /********************************************************************************************/
  /**
   * A standard copy constructor
   *
   * @param cp A copy of another GPsi2STo2K2PiGamIndividual
   */
  GPsi2STo2K2PiGamIndividual(const GPsi2STo2K2PiGamIndividual& cp)
    : GParameterSet(cp)
    ,_psi2STo2K2PiGamLhPtr( cp.psi2STo2K2PiGamLhPtr()->clone_() )
    ,_upar( cp._upar)
  {
    setUp();
  }
  
  /********************************************************************************************/
  /**
   * The standard destructor
   */
  virtual ~GPsi2STo2K2PiGamIndividual()
  { 
  }
  
  /********************************************************************************************/
  
  boost::shared_ptr<AbsPsi2STo2K2PiGamLh> psi2STo2K2PiGamLhPtr() const {return _psi2STo2K2PiGamLhPtr;}

  /**
   * A standard assignment operator
   *
   * @param cp A copy of another GPsi2STo2K2PiGamIndividual object
   * @return A constant reference to this object
   */
  const GPsi2STo2K2PiGamIndividual& operator=(const GPsi2STo2K2PiGamIndividual& cp){
    GPsi2STo2K2PiGamIndividual::load_(&cp);
    return *this;
  }
  
  /*******************************************************************************************/
  /**
   * Checks for equality with another GPsi2STo2K2PiGamIndividual object.
   *
   * NOTE: THIS FUNCTION IS OPTIONAL AND IS MAINLY USED IN CONJUNCTION WITH UNIT TESTS.
   * You do not need it if you do not intend to perform unit tests.
   *
   * @param  cp A constant reference to another GPsi2STo2K2PiGamIndividual object
   * @return A boolean indicating whether both objects are equal
   */
  bool operator==(const GPsi2STo2K2PiGamIndividual& cp) const {
    using namespace Gem::Util;
    // Means: The expectation of equality was fulfilled, if no error text was emitted (which converts to "true")
    return !checkRelationshipWith(cp, CE_EQUALITY, 0.,"GPsi2STo2K2PiGamIndividual::operator==","cp", CE_SILENT);
  }
  
  /*******************************************************************************************/
  /**
   * Checks for inequality with another GPsi2STo2K2PiGamIndividual object.
   *
   * NOTE: THIS FUNCTION IS OPTIONAL AND IS MAINLY USED IN CONJUNCTION WITH UNIT TESTS.
   * You do not need it if you do not intend to perform unit tests.
   *
   * @param  cp A constant reference to another GPsi2STo2K2PiGamIndividual object
   * @return A boolean indicating whether both objects are inequal
   */
  bool operator!=(const GPsi2STo2K2PiGamIndividual& cp) const {
    using namespace Gem::Util;
    // Means: The expectation of inequality was fulfilled, if no error text was emitted (which converts to "true")
    return !checkRelationshipWith(cp, CE_INEQUALITY, 0.,"GPsi2STo2K2PiGamIndividual::operator!=","cp", CE_SILENT);
  }


  /*******************************************************************************************/
  /**
   * Checks whether a given expectation for the relationship between this object and another object
   * is fulfilled.
   *
   * NOTE: THIS FUNCTION IS OPTIONAL AND IS MAINLY USED IN CONJUNCTION WITH UNIT TESTS.
   * You do not need it if you do not intend to perform unit tests.
   *
   * @param cp A constant reference to another object, camouflaged as a GObject
   * @param e The expected outcome of the comparison
   * @param limit The maximum deviation for floating point values (important for similarity checks)
   * @param caller An identifier for the calling entity
   * @param y_name An identifier for the object that should be compared to this one
   * @param withMessages Whether or not information should be emitted in case of deviations from the expected outcome
   * @return A boost::optional<std::string> object that holds a descriptive string if expectations were not met
   */
  boost::optional<std::string> checkRelationshipWith(const GObject& cp,
						     const Gem::Util::expectation& e,
						     const double& limit,
						     const std::string& caller,
						     const std::string& y_name,
						     const bool& withMessages) const
  {
    using namespace Gem::Util;
    using namespace Gem::Util::POD;
    
    // Check that we are not accidently assigning this object to itself
    selfAssignmentCheck<GPsi2STo2K2PiGamIndividual>(&cp);
    // Use this call instead when local data needs to be loaded:
    // const GPsi2STo2K2PiGamIndividual *p_load = GObject::conversion_cast<GPsi2STo2K2PiGamIndividual>(&cp);
    
    // Will hold possible deviations from the expectation, including explanations
    std::vector<boost::optional<std::string> > deviations;
    
    // Check our parent class'es data ...
    deviations.push_back(GParameterSet::checkRelationshipWith(cp, e, limit, "GPsi2STo2K2PiGamIndividual", y_name, withMessages));
    
    // Check local data like this:
    // deviations.push_back(checkExpectation(withMessages, "GPsi2STo2K2PiGamIndividual", val_, p_load->val_, "val_", "p_load->val_", e , limit));
    
    return evaluateDiscrepancies("GPsi2STo2K2PiGamIndividual", caller, deviations, e);
  }
  
  /*******************************************************************************************/
  /**
   * Applies modifications to this object. This is needed for testing purposes
   *
   * @return A boolean which indicates whether modifications were made
   */
  virtual bool modify_GUnitTests() {
    bool result;
    
    // Call the parent class'es function
    if(GParameterSet::modify_GUnitTests()) result = true;
    
    // Check that this individual actually contains data to be modified
    if(this->size() != 0) {
      mutate(); // Perform modifications
      result = true;
    }
    
    return result;
  }
  
  /*******************************************************************************************/
  /**
   * Performs self tests that are expected to succeed. This is needed for testing purposes
   */
  virtual void specificTestsNoFailureExpected_GUnitTests() {
    const boost::uint32_t NITERATIONS=100;
    
    // Call the parent class'es function
    GParameterSet::specificTestsNoFailureExpected_GUnitTests();
    
    // Create an individual
    boost::shared_ptr<Gem::GenEvA::GPsi2STo2K2PiGamIndividual> p
      = boost::shared_ptr<Gem::GenEvA::GPsi2STo2K2PiGamIndividual>(new GPsi2STo2K2PiGamIndividual(*this));
    
    // Mutate a number of times and check that there were changes
    double oldfitness = p->fitness();
    for(boost::uint32_t i=0; i<NITERATIONS; i++) {
      p->mutate();
      double newfitness = p->fitness();
      BOOST_CHECK_MESSAGE(newfitness != oldfitness, "Rare failures are normal for this test / " << i << "/" << NITERATIONS);
      oldfitness = newfitness;
    }
  }
  
  /*******************************************************************************************/
  /**
   * Performs self tests that are expected to fail. This is needed for testing purposes
   */
  virtual void specificTestsFailuresExpected_GUnitTests() {
    // Call the parent class'es function
    GParameterSet::specificTestsFailuresExpected_GUnitTests();
  }

  bool getFitParams(Psi2STo2K2PiGamData::fitParamVal& fitParmVal)
  { 
    boost::shared_ptr<GBoundedDoubleCollection> vC = pc_at<GBoundedDoubleCollection>(0);

    std::map<size_t, double>::const_iterator itFixedParams;    
    std::vector<double> par;
    for(std::size_t i=0; i<vC->size(); i++){

      if ( !(_fixedParams.find(i) == _fixedParams.end()) ){
	par.push_back(_fixedParams[i]);
      }

      double value = vC->at(i)->value();
      par.push_back(value);
    }



    _psi2STo2K2PiGamLhPtr->setFitParamVal(fitParmVal, par);
    
    return true;
  }
  
  void printFitParams(Psi2STo2K2PiGamData::fitParamVal& theParamVal){
    _psi2STo2K2PiGamLhPtr->printCurrentFitResult(theParamVal);
  }  

protected:
  /********************************************************************************************/
  /**
   * Loads the data of another GPsi2STo2K2PiGamIndividual, camouflaged as a GObject.
   *
   * @param cp A copy of another GPsi2STo2K2PiGamIndividual, camouflaged as a GObject
   */
  virtual void load_(const GObject* cp)
  {
    // Check that we are not accidently assigning this object to itself
    selfAssignmentCheck<GPsi2STo2K2PiGamIndividual>(cp);
    // Use this call instead when local data needs to be loaded:
    // const GPsi2STo2K2PiGamIndividual *p_load = GObject::conversion_cast<GPsi2STo2K2PiGamIndividual>(cp);
    
    // Load our parent's data
    GParameterSet::load_(cp);
    
    // Load local data here like this:
    // myVar = p_load->myVar;
  }
  
  /********************************************************************************************/
  /**
   * Creates a deep clone of this object
   *
   * @return A deep clone of this object, camouflaged as a GObject
   */
  virtual GObject* clone_() const {
    return new GPsi2STo2K2PiGamIndividual(*this);
  }
  
  /********************************************************************************************/
  /**
   * The actual fitness calculation takes place here.
   *
   * @return The value of this object
   */
  virtual double fitnessCalculation(){
    double result = 0.;

    Psi2STo2K2PiGamData::fitParamVal currentFitParms;
    getFitParams(currentFitParms);
    result=_psi2STo2K2PiGamLhPtr->calcLogLh(currentFitParms);

    return result;
  }


private:
  /********************************************************************************************/
  /**
   * The default constructor. Intentionally private and empty, as it is only needed for
   * serialization purposes.
   */
  boost::shared_ptr<AbsPsi2STo2K2PiGamLh> _psi2STo2K2PiGamLhPtr;
  MnUserParameters _upar;
  std::map<size_t, double> _fixedParams;
 

  GPsi2STo2K2PiGamIndividual() :GParameterSet()
  {	/* nothing */ }


  void setUp(){
		// Set up a GBoundedDoubleCollection
		boost::shared_ptr<GBoundedDoubleCollection> gbdc_ptr(new GBoundedDoubleCollection());

		// Create a suitable adaptor (sigma=0.1, sigma-adaption=0.5, min sigma=0, max sigma=0,5)
		boost::shared_ptr<GDoubleGaussAdaptor> gdga_ptr(new GDoubleGaussAdaptor(0.1, 0.5, 0., 0.5));
		gdga_ptr->setAdaptionThreshold(1); // Mutation parameters are adapted after each mutation
		gdga_ptr->setRnrGenerationMode(Gem::Util::RNRFACTORY); // Random number generation in the factory
// 		gdga_ptr->setMutationProbability(0.05); // The likelihood for a parameter to be mutated
		gdga_ptr->setMutationProbability(0.2); // The likelihood for a parameter to be mutated

		// Register the adaptor with the collection. You could also add individual adaptors
		// to the GBoundedDouble objects below.
		gbdc_ptr->addAdaptor(gdga_ptr);
		
		const std::vector<ROOT::Minuit2::MinuitParameter> theMinuitParameters= _upar.Parameters();

		size_t tmpCounter=0;		
		std::vector<ROOT::Minuit2::MinuitParameter>::const_iterator it;
		for ( it=theMinuitParameters.begin(); it != theMinuitParameters.end(); ++it){
		  double lowerLimit=-100000.;
		  double upperLimit=100000.;

		  if (it->IsFixed()){
		    _fixedParams[tmpCounter]=it->Value();
                    tmpCounter++;
		    continue;
		  }

		  if (it->HasLowerLimit()) lowerLimit=it->LowerLimit();
		  if (it->HasUpperLimit()) upperLimit=it->UpperLimit();
		  
		  boost::shared_ptr<GBoundedDouble> gbd_ptr(new GBoundedDouble(lowerLimit, upperLimit) );
		  gbdc_ptr->push_back(gbd_ptr);
		  tmpCounter++; 
		}
		
		// Add the collection to this object
		this->push_back(gbdc_ptr);

  }

  /********************************************************************************************/
  // You can add other variables here. Do not forget to serialize them if necessary
  // int myVar;
};
  
  /*************************************************************************************************/
  
} /* namespace GenEvA */
} /* namespace Gem */

// Needed for serialization purposes
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT(Gem::GenEvA::GPsi2STo2K2PiGamIndividual)
  
// Needed for testing purposes
/*************************************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////
/*************************************************************************************************/
/**
 * As the Gem::GenEvA::Gem::GenEvA::GPsi2STo2K2PiGamIndividual has a private default constructor, we need to provide a
 * specialization of the factory function that creates GStartProjectIndividual objects
 */
// template <>
// boost::shared_ptr<Gem::GenEvA::GPsi2STo2K2PiGamIndividual> TFactory_GUnitTests<Gem::GenEvA::GPsi2STo2K2PiGamIndividual>() {
// 	return boost::shared_ptr<Gem::GenEvA::GPsi2STo2K2PiGamIndividual>(new Gem::GenEvA::GPsi2STo2K2PiGamIndividual(1000,-10.,10.));
// }

/*************************************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////
/*************************************************************************************************/

#endif /* GSTARTINDIVIDUAL_HPP_ */

