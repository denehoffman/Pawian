// Standard header files go here
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm> // for std::sort
#include <utility> // For std::pair

// Includes check for correct Boost version(s)
#include "GGlobalDefines.hpp"

// Boost header files go here
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/cast.hpp>

#ifndef GEtacToapi0Individual_HPP_
#define GEtacToapi0Individual_HPP_

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

#include "ErrLogger/ErrLineLog.hh"

#include "Examples/EtacToapi0Fit/EtacToapi0EventList.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
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
class GEtacToapi0Individual :public GParameterSet
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
	GEtacToapi0Individual(boost::shared_ptr<const EtacToapi0EventList> theEvtList)
	  : GParameterSet()
	  ,_etacLhPtr( new EtacToapi0Lh(theEvtList) )
	  {
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

		// Add bounded double objects
		boost::shared_ptr<GBoundedDouble> gbd_ptr(new GBoundedDouble(0.7, 2.3) );  //a mass
		gbdc_ptr->push_back(gbd_ptr);
		boost::shared_ptr<GBoundedDouble> gbd_ptr1(new GBoundedDouble(0.01, 0.8) ); //a width
		gbdc_ptr->push_back(gbd_ptr1);
		boost::shared_ptr<GBoundedDouble> gbd_ptr2(new GBoundedDouble(0., 1.) ); //spin0 content
		gbdc_ptr->push_back(gbd_ptr2);
		boost::shared_ptr<GBoundedDouble> gbd_ptr3(new GBoundedDouble(0., 1.) ); //spin1 content
		gbdc_ptr->push_back(gbd_ptr3);
		boost::shared_ptr<GBoundedDouble> gbd_ptr4(new GBoundedDouble(0., 1.) ); //spin2 content
		gbdc_ptr->push_back(gbd_ptr4);
		
		
		// Add the collection to this object
		this->push_back(gbdc_ptr);
	  }
  
  /********************************************************************************************/
  /**
   * A standard copy constructor
   *
   * @param cp A copy of another GEtacToapi0Individual
   */
  GEtacToapi0Individual(const GEtacToapi0Individual& cp)
    : GParameterSet(cp)
    ,_etacLhPtr( new EtacToapi0Lh(cp.etacLhPtr()) )
  { /* nothing */ }
  
  /********************************************************************************************/
  /**
   * The standard destructor
   */
  virtual ~GEtacToapi0Individual()
  { 
  }
  
  /********************************************************************************************/
  
  boost::shared_ptr<EtacToapi0Lh> etacLhPtr() const {return _etacLhPtr;}
  /**
   * A standard assignment operator
   *
   * @param cp A copy of another GEtacToapi0Individual object
   * @return A constant reference to this object
   */
  const GEtacToapi0Individual& operator=(const GEtacToapi0Individual& cp){
    GEtacToapi0Individual::load_(&cp);
    return *this;
  }
  
  /*******************************************************************************************/
  /**
   * Checks for equality with another GEtacToapi0Individual object.
   *
   * NOTE: THIS FUNCTION IS OPTIONAL AND IS MAINLY USED IN CONJUNCTION WITH UNIT TESTS.
   * You do not need it if you do not intend to perform unit tests.
   *
   * @param  cp A constant reference to another GEtacToapi0Individual object
   * @return A boolean indicating whether both objects are equal
   */
  bool operator==(const GEtacToapi0Individual& cp) const {
    using namespace Gem::Util;
    // Means: The expectation of equality was fulfilled, if no error text was emitted (which converts to "true")
    return !checkRelationshipWith(cp, CE_EQUALITY, 0.,"GEtacToapi0Individual::operator==","cp", CE_SILENT);
  }
  
  /*******************************************************************************************/
  /**
   * Checks for inequality with another GEtacToapi0Individual object.
   *
   * NOTE: THIS FUNCTION IS OPTIONAL AND IS MAINLY USED IN CONJUNCTION WITH UNIT TESTS.
   * You do not need it if you do not intend to perform unit tests.
   *
   * @param  cp A constant reference to another GEtacToapi0Individual object
   * @return A boolean indicating whether both objects are inequal
   */
  bool operator!=(const GEtacToapi0Individual& cp) const {
    using namespace Gem::Util;
    // Means: The expectation of inequality was fulfilled, if no error text was emitted (which converts to "true")
    return !checkRelationshipWith(cp, CE_INEQUALITY, 0.,"GEtacToapi0Individual::operator!=","cp", CE_SILENT);
  }


  bool getFitParams(fitParamVal& fitParmVal)
  { 
    std::vector<double> theParms;
    // Extract the GDoubleCollection object. In a realistic scenario, you might want
    // to add error checks here upon first invocation.
    boost::shared_ptr<GBoundedDoubleCollection> vC = pc_at<GBoundedDoubleCollection>(0);

    for(std::size_t i=0; i<vC->size(); i++){
      double value = vC->at(i)->value();
      theParms.push_back(value);
    }
    
    return _etacLhPtr->setFitParamVal(fitParmVal, theParms);
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
    selfAssignmentCheck<GEtacToapi0Individual>(&cp);
    // Use this call instead when local data needs to be loaded:
    // const GEtacToapi0Individual *p_load = GObject::conversion_cast<GEtacToapi0Individual>(&cp);
    
    // Will hold possible deviations from the expectation, including explanations
    std::vector<boost::optional<std::string> > deviations;
    
    // Check our parent class'es data ...
    deviations.push_back(GParameterSet::checkRelationshipWith(cp, e, limit, "GEtacToapi0Individual", y_name, withMessages));
    
    // Check local data like this:
    // deviations.push_back(checkExpectation(withMessages, "GEtacToapi0Individual", val_, p_load->val_, "val_", "p_load->val_", e , limit));
    
    return evaluateDiscrepancies("GEtacToapi0Individual", caller, deviations, e);
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
    boost::shared_ptr<Gem::GenEvA::GEtacToapi0Individual> p
      = boost::shared_ptr<Gem::GenEvA::GEtacToapi0Individual>(new GEtacToapi0Individual(*this));
    
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
  
protected:
  /********************************************************************************************/
  /**
   * Loads the data of another GEtacToapi0Individual, camouflaged as a GObject.
   *
   * @param cp A copy of another GEtacToapi0Individual, camouflaged as a GObject
   */
  virtual void load_(const GObject* cp)
  {
    // Check that we are not accidently assigning this object to itself
    selfAssignmentCheck<GEtacToapi0Individual>(cp);
    // Use this call instead when local data needs to be loaded:
    // const GEtacToapi0Individual *p_load = GObject::conversion_cast<GEtacToapi0Individual>(cp);
    
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
    return new GEtacToapi0Individual(*this);
  }
  
  /********************************************************************************************/
  /**
   * The actual fitness calculation takes place here.
   *
   * @return The value of this object
   */
  virtual double fitnessCalculation(){
    double result = 0.;
    
    fitParamVal theFitParmValTmp;
    assert(getFitParams(theFitParmValTmp));   
    result=_etacLhPtr->calcLogLh(theFitParmValTmp);
    ErrMsg(debugging) << "current LH= " << result << endmsg;    
    return result;
  }
  
private:
  /********************************************************************************************/
  /**
   * The default constructor. Intentionally private and empty, as it is only needed for
   * serialization purposes.
   */
  boost::shared_ptr<EtacToapi0Lh> _etacLhPtr;

  GEtacToapi0Individual() :GParameterSet()
  {	/* nothing */ }
  
  /********************************************************************************************/
  // You can add other variables here. Do not forget to serialize them if necessary
  // int myVar;
};
  
  /*************************************************************************************************/
  
} /* namespace GenEvA */
} /* namespace Gem */

// Needed for serialization purposes
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT(Gem::GenEvA::GEtacToapi0Individual)
  
// Needed for testing purposes
/*************************************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////
/*************************************************************************************************/
/**
 * As the Gem::GenEvA::Gem::GenEvA::GEtacToapi0Individual has a private default constructor, we need to provide a
 * specialization of the factory function that creates GStartProjectIndividual objects
 */
// template <>
// boost::shared_ptr<Gem::GenEvA::GEtacToapi0Individual> TFactory_GUnitTests<Gem::GenEvA::GEtacToapi0Individual>() {
// 	return boost::shared_ptr<Gem::GenEvA::GEtacToapi0Individual>(new Gem::GenEvA::GEtacToapi0Individual(1000,-10.,10.));
// }

/*************************************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////
/*************************************************************************************************/

#endif /* GSTARTINDIVIDUAL_HPP_ */

