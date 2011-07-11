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

#ifndef GOMEGAPIINDIVIDUALLS_HPP_
#define GOMEGAPIINDIVIDUALLS_HPP_

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

#include "Examples/pbarpToOmegaPiLS/OmegaPiDataLS.hh"
// #include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
// #include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiLhLS.hh"

namespace Gem
{
  namespace Pawian
  {
    /************************************************************************************************/
    /**
     * This individual performs a partial wave analysis using the Geneva library
     */
    class GOmegaPiIndividualLS 
      :public Gem::Geneva::GParameterSet
    {
      ///////////////////////////////////////////////////////////////////////
      friend class boost::serialization::access;

      template<typename Archive>
      void serialize(Archive & ar, const unsigned int) {
	using boost::serialization::make_nvp;

	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Gem::Geneva::GParameterSet);
      }
      ///////////////////////////////////////////////////////////////////////

    public:
      /** @brief A constructor that initializes the object with a likelihood calculator */
      GOmegaPiIndividualLS(boost::shared_ptr<AbsOmegaPiLhLS>);
      /** @brief The standard copy constructor */
      GOmegaPiIndividualLS(const GOmegaPiIndividualLS&);
      /** @brief The destructor */
      virtual ~GOmegaPiIndividualLS();

      /** @brief Loads user-specified data */
      virtual void loadConstantData(boost::shared_ptr<GObject>);
      
      /** @brief A standard assignment operator */
      const GOmegaPiIndividualLS& operator=(const GOmegaPiIndividualLS&);

      /** @brief Gives access to the likelihood calculator */
      boost::shared_ptr<AbsOmegaPiLhLS> omegaPiLhPtr() const;
      /** @brief Gives access to the fit parameters */
      bool getFitParams(OmegaPiDataLS::fitParamVal&);
      /** @brief Prints the fit parameters */
      void printFitParams(OmegaPiDataLS::fitParamVal&);
  
    protected:
      /** @brief Loads the data of another GOmegaPiIndividualLS */
      virtual void load_(const Gem::Geneva::GObject*);
      /** Creates a deep clone of this object */
      virtual Gem::Geneva::GObject* clone_() const;

      /** @brief The actual fitness calculation takes place here */
      virtual double fitnessCalculation();
  
    private:
      /** @brief The default constructor. Intentionally private */
      GOmegaPiIndividualLS();

       
      boost::shared_ptr<AbsOmegaPiLhLS> _omegaPiLhPtr; ///< Holds the likelihood calculator
    };
  
    /*************************************************************************************************/
  
  } /* namespace Pawian */
} /* namespace Gem */

BOOST_CLASS_EXPORT_KEY(Gem::Pawian::GOmegaPiIndividualLS)

#endif /* GOMEGAPIINDIVIDUALLS_HPP_ */
