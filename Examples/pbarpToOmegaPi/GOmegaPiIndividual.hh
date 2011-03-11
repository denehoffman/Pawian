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
#include <common/GGlobalDefines.hpp>

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
#include <hap/GRandomT.hpp>
#include <common/GCommonEnums.hpp>
#include <common/GExceptions.hpp>
#include <geneva/GConstrainedDoubleObject.hpp>
#include <geneva/GConstrainedDoubleObjectCollection.hpp>
#include <geneva/GDoubleGaussAdaptor.hpp>
#include <geneva/GObjectExpectationChecksT.hpp>
#include <geneva/GParameterObjectCollection.hpp>
#include <geneva/GParameterSet.hpp>

#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"

#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "Examples/pbarpToOmegaPi/AbsOmegaPiLh.hh"

namespace Gem
{
  namespace Geneva
  {
    /************************************************************************************************/
    /**
     * This individual performs a partial wave analysis using the Geneva library
     */
    class GOmegaPiIndividual 
      :public GParameterSet
    {
      ///////////////////////////////////////////////////////////////////////
      friend class boost::serialization::access;

      template<typename Archive>
      void serialize(Archive & ar, const unsigned int) {
	using boost::serialization::make_nvp;

	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GParameterSet);
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
	// Set up a GConstrainedDoubleObjectCollection
	boost::shared_ptr<GConstrainedDoubleObjectCollection> gbdc_ptr(new GConstrainedDoubleObjectCollection());

	// Create a suitable adaptor (sigma=0.1, sigma-adaption=0.5, min sigma=0, max sigma=0,5)
	boost::shared_ptr<GDoubleGaussAdaptor> gdga_ptr(new GDoubleGaussAdaptor(0.1, 0.5, 0., 0.5));

	gdga_ptr->setAdaptionThreshold(1); // Adaption parameters are modified after each adaption
	gdga_ptr->setAdaptionProbability(0.05); // The likelihood for a parameter to be adapted

	boost::shared_ptr<const pbarpToOmegaPi0States> theStates=_omegaPiLhPtr->omegaPi0States();

	std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=theStates->jpclsSinglet();
	setFitParamVal(JPCLSOmegaSinglet, gbdc_ptr);

	std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=theStates->jpclsTriplet0();
	setFitParamVal(JPCLSOmegaTriplet0, gbdc_ptr);

	std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=theStates->jpclsTriplet1();
	setFitParamVal(JPCLSOmegaTriplet1, gbdc_ptr);


	for(std::size_t i=0; i<gbdc_ptr->size(); i++) {
	  // Register the adaptor with GConstrainedDoubleObject objects
	  gbdc_ptr->at(i)->addAdaptor(gdga_ptr);
	}

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
      { /* nothing */ }
  
      /********************************************************************************************/
      /**
       * The standard destructor
       */
      virtual ~GOmegaPiIndividual()
      { /* nothing */ }
  
      /********************************************************************************************/
  
      boost::shared_ptr<AbsOmegaPiLh> omegaPiLhPtr() const {
	return _omegaPiLhPtr;
      }

      /********************************************************************************************/
      /**
       * A standard assignment operator
       *
       * @param cp A copy of another GOmegaPiIndividual object
       * @return A constant reference to this object
       */
      const GOmegaPiIndividual& operator=(const GOmegaPiIndividual& cp) {
	GOmegaPiIndividual::load_(&cp);
	return *this;
      }
  
      /********************************************************************************************/

      bool getFitParams(OmegaPiData::fitParamVal& fitParmVal) { 
	std::vector<double> theParms;
	// Extract the GDoubleCollection object. In a realistic scenario, you might want
	// to add error checks here upon first invocation.
	boost::shared_ptr<GConstrainedDoubleObjectCollection> vC = at<GConstrainedDoubleObjectCollection>(0);
    
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
	for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS) {
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
    
	for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS) {
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
    
	for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS) {
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

      /********************************************************************************************/

      void printFitParams(OmegaPiData::fitParamVal& fitParmVal) {

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

      /********************************************************************************************/

      boost::shared_ptr<AbsOmegaPiLh> getOmegaPiLhPtr() {
	return _omegaPiLhPtr;
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
      GOmegaPiIndividual() :GParameterSet()
      {	/* nothing */ }


      /********************************************************************************************/
      void setFitParamVal(std::vector< boost::shared_ptr<const JPCLS> > theJPCLSs, boost::shared_ptr<GConstrainedDoubleObjectCollection> theGbdc_ptr){
    
	std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
	int counter=0;
	for ( itJPCLS=theJPCLSs.begin(); itJPCLS!=theJPCLSs.end(); ++itJPCLS){
	  //now fill the fitParameterMap
      
	  boost::shared_ptr<GConstrainedDoubleObject> gbd_ptr(new GConstrainedDoubleObject(0., 1.) ); //JPCLS magnitude
	  theGbdc_ptr->push_back(gbd_ptr);

	  if (counter>0){ 
	    boost::shared_ptr<GConstrainedDoubleObject>  gbd_ptr(new GConstrainedDoubleObject(-M_PI, M_PI) ); //JPCLS phi
	    theGbdc_ptr->push_back(gbd_ptr);
	  }
	  counter++; 
	}
      }

      /********************************************************************************************/
      boost::shared_ptr<AbsOmegaPiLh> _omegaPiLhPtr;
    };
  
    /*************************************************************************************************/
  
  } /* namespace Geneva */
} /* namespace Gem */

// Needed for serialization purposes
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT(Gem::Geneva::GOmegaPiIndividual)

#endif /* GSTARTINDIVIDUAL_HPP_ */
