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

#ifndef GOmegaPiIndividual_HPP_
#define GOmegaPiIndividual_HPP_

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

#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"

#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "Examples/pbarpToOmegaPi/AbsOmegaPiLh.hh"

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
class GOmegaPiIndividual :public GParameterSet
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
	GOmegaPiIndividual(boost::shared_ptr<AbsOmegaPiLh> theLh)
	  : GParameterSet()
	  ,_omegaPiLhPtr( theLh->clone_() )
	  {
		// Set up a GBoundedDoubleCollection
		boost::shared_ptr<GBoundedDoubleCollection> gbdc_ptr(new GBoundedDoubleCollection());

		// Create a suitable adaptor (sigma=0.1, sigma-adaption=0.5, min sigma=0, max sigma=0,5)
		boost::shared_ptr<GDoubleGaussAdaptor> gdga_ptr(new GDoubleGaussAdaptor(0.1, 0.5, 0., 0.5));
		gdga_ptr->setAdaptionThreshold(1); // Mutation parameters are adapted after each mutation
		gdga_ptr->setRnrGenerationMode(Gem::Util::RNRFACTORY); // Random number generation in the factory
// 		gdga_ptr->setMutationProbability(0.05); // The likelihood for a parameter to be mutated
		gdga_ptr->setMutationProbability(0.1); // The likelihood for a parameter to be mutated

		// Register the adaptor with the collection. You could also add individual adaptors
		// to the GBoundedDouble objects below.
		gbdc_ptr->addAdaptor(gdga_ptr);

                
		boost::shared_ptr<const pbarpToOmegaPi0States> theStates=_omegaPiLhPtr->omegaPi0States();

		std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=theStates->jpclsSinglet();
		setFitParamVal(JPCLSOmegaSinglet, gbdc_ptr);

		std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=theStates->jpclsTriplet0();
		setFitParamVal(JPCLSOmegaTriplet0, gbdc_ptr);

		std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=theStates->jpclsTriplet1();
		setFitParamVal(JPCLSOmegaTriplet1, gbdc_ptr);

		// Add the collection to this object
		this->push_back(gbdc_ptr);
	  }
  
  /********************************************************************************************/
  /**
   * A standard copy constructor
   *
   * @param cp A copy of another GOmegaPiIndividual
   */
  GOmegaPiIndividual(const GOmegaPiIndividual& cp)
    : GParameterSet(cp)
    ,_omegaPiLhPtr( cp.omegaPiLhPtr()->clone_() )
  { 
  }
  
  /********************************************************************************************/
  /**
   * The standard destructor
   */
  virtual ~GOmegaPiIndividual()
  { 
  }
  
  /********************************************************************************************/
  
  boost::shared_ptr<AbsOmegaPiLh> omegaPiLhPtr() const {return _omegaPiLhPtr;}
  /**
   * A standard assignment operator
   *
   * @param cp A copy of another GOmegaPiIndividual object
   * @return A constant reference to this object
   */
  const GOmegaPiIndividual& operator=(const GOmegaPiIndividual& cp){
    GOmegaPiIndividual::load_(&cp);
    return *this;
  }
  
  /*******************************************************************************************/
  /**
   * Checks for equality with another GOmegaPiIndividual object.
   *
   * NOTE: THIS FUNCTION IS OPTIONAL AND IS MAINLY USED IN CONJUNCTION WITH UNIT TESTS.
   * You do not need it if you do not intend to perform unit tests.
   *
   * @param  cp A constant reference to another GOmegaPiIndividual object
   * @return A boolean indicating whether both objects are equal
   */
  bool operator==(const GOmegaPiIndividual& cp) const {
    using namespace Gem::Util;
    // Means: The expectation of equality was fulfilled, if no error text was emitted (which converts to "true")
    return !checkRelationshipWith(cp, CE_EQUALITY, 0.,"GOmegaPiIndividual::operator==","cp", CE_SILENT);
  }
  
  /*******************************************************************************************/
  /**
   * Checks for inequality with another GOmegaPiIndividual object.
   *
   * NOTE: THIS FUNCTION IS OPTIONAL AND IS MAINLY USED IN CONJUNCTION WITH UNIT TESTS.
   * You do not need it if you do not intend to perform unit tests.
   *
   * @param  cp A constant reference to another GOmegaPiIndividual object
   * @return A boolean indicating whether both objects are inequal
   */
  bool operator!=(const GOmegaPiIndividual& cp) const {
    using namespace Gem::Util;
    // Means: The expectation of inequality was fulfilled, if no error text was emitted (which converts to "true")
    return !checkRelationshipWith(cp, CE_INEQUALITY, 0.,"GOmegaPiIndividual::operator!=","cp", CE_SILENT);
  }


  bool getFitParams(OmegaPiData::fitParamVal& fitParmVal)
  { 
    std::vector<double> theParms;
    // Extract the GDoubleCollection object. In a realistic scenario, you might want
    // to add error checks here upon first invocation.
    boost::shared_ptr<GBoundedDoubleCollection> vC = pc_at<GBoundedDoubleCollection>(0);
    
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

    boost::shared_ptr<const pbarpToOmegaPi0States> theStates=_omegaPiLhPtr->omegaPi0States();

    std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=theStates->jpclsSinglet();
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=theStates->jpclsTriplet0();
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=theStates->jpclsTriplet1();
    
    
    std::vector<double> par;
    for(std::size_t i=0; i<vC->size(); i++){
      double value = vC->at(i)->value();
      par.push_back(value);
    }
    
    if (par.size()!= JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2+JPCLSOmegaTriplet1.size()*2-3){
      std::cout << "size of parameters wrong!!! par.size()=" << par.size() << "\tJPCLSOmegaSinglet.size()+JPCLSOmegaTriplet0.size()+JPCLSOmegaTriplet1.size()-3=" << JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2+JPCLSOmegaTriplet1.size()*2-3 << std::endl;
      assert(0);
    }

    
    unsigned int counter=0;
    for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
      //now fill the fitParameterMap
      double mag=par[counter];
      counter++;
      double phi=0.;
      if (counter>1){ phi=par[counter];
	counter++;
      }
      std::pair <double,double> tmpParameter=make_pair(mag,phi);
      fitParmVal.omegaProdSinglet[(*itJPCLS)]=tmpParameter; 
    }
    
    for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
      //now fill the fitParameterMap
      double mag=par[counter];
      counter++;
      double phi=0.;
      if ( counter > JPCLSOmegaSinglet.size()*2){ phi=par[counter];
	counter++;
      }
      std::pair <double,double> tmpParameter=make_pair(mag,phi);
      fitParmVal.omegaProdTriplet0[(*itJPCLS)]=tmpParameter; 
    }
    
    for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
      //now fill the fitParameterMap
      double mag=par[counter];
      counter++;
      double phi=0.;
      if (counter>JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2-1){ phi=par[counter];
	counter++;
      }
      std::pair <double,double> tmpParameter=make_pair(mag,phi);
      fitParmVal.omegaProdTriplet1[(*itJPCLS)]=tmpParameter; 
    }
    
    return true;
  }

 void printFitParams(OmegaPiData::fitParamVal& fitParmVal)
  {

    boost::shared_ptr<const pbarpToOmegaPi0States> theStates=_omegaPiLhPtr->omegaPi0States();

    std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=theStates->jpclsSinglet();
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=theStates->jpclsTriplet0();
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=theStates->jpclsTriplet1();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

    std::cout << "***fit parameter singlet states*** " <<std::endl;  
    for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
      std::cout << (*itJPCLS)->name()<< "\t";
      std::pair<double, double> tmpParam=fitParmVal.omegaProdSinglet[(*itJPCLS)];
      std::cout <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << std::endl;
    }

    std::cout << "***fit parameter triplet m=0 states*** " <<std::endl;  
    for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
      std::cout<< (*itJPCLS)->name()<< "\t";
      std::pair<double, double> tmpParam=fitParmVal.omegaProdTriplet0[(*itJPCLS)];
      std::cout <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << std::endl;
    }
    
    std::cout << "***fit parameter triplet m=1 states*** " <<std::endl;  
    for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
      std::cout<< (*itJPCLS)->name()<< "\t";
      std::pair<double, double> tmpParam=fitParmVal.omegaProdTriplet1[(*itJPCLS)];
      std::cout <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << std::endl;
    }
 
    std::cout << std::endl;
    return;
  }

  boost::shared_ptr<AbsOmegaPiLh> getOmegaPiLhPtr() {return _omegaPiLhPtr;}
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
    selfAssignmentCheck<GOmegaPiIndividual>(&cp);
    // Use this call instead when local data needs to be loaded:
    // const GOmegaPiIndividual *p_load = GObject::conversion_cast<GOmegaPiIndividual>(&cp);
    
    // Will hold possible deviations from the expectation, including explanations
    std::vector<boost::optional<std::string> > deviations;
    
    // Check our parent class'es data ...
    deviations.push_back(GParameterSet::checkRelationshipWith(cp, e, limit, "GOmegaPiIndividual", y_name, withMessages));
    
    // Check local data like this:
    // deviations.push_back(checkExpectation(withMessages, "GOmegaPiIndividual", val_, p_load->val_, "val_", "p_load->val_", e , limit));
    
    return evaluateDiscrepancies("GOmegaPiIndividual", caller, deviations, e);
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
    boost::shared_ptr<Gem::GenEvA::GOmegaPiIndividual> p
      = boost::shared_ptr<Gem::GenEvA::GOmegaPiIndividual>(new GOmegaPiIndividual(*this));
    
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
   * Loads the data of another GOmegaPiIndividual, camouflaged as a GObject.
   *
   * @param cp A copy of another GOmegaPiIndividual, camouflaged as a GObject
   */
  virtual void load_(const GObject* cp)
  {
    // Check that we are not accidently assigning this object to itself
    selfAssignmentCheck<GOmegaPiIndividual>(cp);
    // Use this call instead when local data needs to be loaded:
    // const GOmegaPiIndividual *p_load = GObject::conversion_cast<GOmegaPiIndividual>(cp);
    
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
    return new GOmegaPiIndividual(*this);
  }
  
  /********************************************************************************************/
  /**
   * The actual fitness calculation takes place here.
   *
   * @return The value of this object
   */
  virtual double fitnessCalculation(){
    double result = 0.;
    
    OmegaPiData::fitParamVal theFitParmValTmp;
    assert(getFitParams(theFitParmValTmp));   
    result=_omegaPiLhPtr->calcLogLh(theFitParmValTmp);

    return result;
  }
  
private:
  /********************************************************************************************/
  /**
   * The default constructor. Intentionally private and empty, as it is only needed for
   * serialization purposes.
   */
  boost::shared_ptr<AbsOmegaPiLh> _omegaPiLhPtr;

  GOmegaPiIndividual() :GParameterSet()
  {	/* nothing */ }



  void setFitParamVal(std::vector< boost::shared_ptr<const JPCLS> > theJPCLSs, boost::shared_ptr<GBoundedDoubleCollection> theGbdc_ptr){
    
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    int counter=0;
    for ( itJPCLS=theJPCLSs.begin(); itJPCLS!=theJPCLSs.end(); ++itJPCLS){
      //now fill the fitParameterMap
      
      boost::shared_ptr<GBoundedDouble> gbd_ptr(new GBoundedDouble(0., 1.) ); //JPCLS magnitude
      theGbdc_ptr->push_back(gbd_ptr);

      if (counter>0){ 
      boost::shared_ptr<GBoundedDouble>  gbd_ptr(new GBoundedDouble(-M_PI, M_PI) ); //JPCLS phi
	theGbdc_ptr->push_back(gbd_ptr);
      }
      counter++; 
    }
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
BOOST_CLASS_EXPORT(Gem::GenEvA::GOmegaPiIndividual)
  
// Needed for testing purposes
/*************************************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////
/*************************************************************************************************/
/**
 * As the Gem::GenEvA::Gem::GenEvA::GOmegaPiIndividual has a private default constructor, we need to provide a
 * specialization of the factory function that creates GStartProjectIndividual objects
 */
// template <>
// boost::shared_ptr<Gem::GenEvA::GOmegaPiIndividual> TFactory_GUnitTests<Gem::GenEvA::GOmegaPiIndividual>() {
// 	return boost::shared_ptr<Gem::GenEvA::GOmegaPiIndividual>(new Gem::GenEvA::GOmegaPiIndividual(1000,-10.,10.));
// }

/*************************************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////
/*************************************************************************************************/

#endif /* GSTARTINDIVIDUAL_HPP_ */

