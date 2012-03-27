// Standard header files go here

#include "Examples/pbarpToOmegaPiLS/GOmegaPiIndividualLS.hh"

BOOST_CLASS_EXPORT_IMPLEMENT(Gem::Pawian::GOmegaPiIndividualLS)

namespace Gem
{
  namespace Pawian
  {
    /********************************************************************************************/
    /**
     * A constructor that initializes the object with a likelihood calculator
     *
     * @param theLh The likelihood calculator
     */
    GOmegaPiIndividualLS::GOmegaPiIndividualLS(boost::shared_ptr<AbsOmegaPiLhLS> theLh)
      : GParameterSet()
      ,_omegaPiLhPtr( theLh->clone_() )
    {
      // Set up a GConstrainedDoubleObjectCollection
      boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> gbdc_ptr(new Gem::Geneva::GConstrainedDoubleObjectCollection());

      _omegaPiLhPtr->setGenevaFitParamVal( gbdc_ptr );

      for(std::size_t i=0; i<gbdc_ptr->size(); i++) {
	// Create a suitable adaptor (sigma=0.1, sigma-adaption=0.5, min sigma=0, max sigma=0,5)
	boost::shared_ptr<Gem::Geneva::GDoubleGaussAdaptor> gdga_ptr(new Gem::Geneva::GDoubleGaussAdaptor(0.1, 0.5, 0., 0.5));
	gdga_ptr->setAdaptionThreshold(1); // Adaption parameters are modified after each adaption
	gdga_ptr->setAdaptionProbability(0.2); // The likelihood for a parameter to be adapted
	
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
     * @param cp A copy of another GOmegaPiIndividualLS
     */
    GOmegaPiIndividualLS::GOmegaPiIndividualLS(const GOmegaPiIndividualLS& cp)
      : GParameterSet(cp)
      , _omegaPiLhPtr( cp._omegaPiLhPtr?boost::shared_ptr<AbsOmegaPiLhLS>((cp._omegaPiLhPtr)->clone_()):boost::shared_ptr<AbsOmegaPiLhLS>() )
    { /* nothing */ }
  
    /********************************************************************************************/
    /**
     * The standard destructor
     */
    GOmegaPiIndividualLS::~GOmegaPiIndividualLS()
    { /* nothing */ }
  
    /********************************************************************************************/
    /** 
     * Loads user-specified data 
     */
    void GOmegaPiIndividualLS::loadConstantData(boost::shared_ptr<GObject> cD_ptr) {
      // Convert GObject smart pointer to local format
      boost::shared_ptr<GOmegaPiIndividualLS> p_load = GObject::conversion_cast<GOmegaPiIndividualLS>(cD_ptr);

      // Load the static data from the data model
      if(p_load->omegaPiLhPtr()) {
	_omegaPiLhPtr = boost::shared_ptr<AbsOmegaPiLhLS>(p_load->omegaPiLhPtr()->clone_());
      }
    }

    /********************************************************************************************/
  
    boost::shared_ptr<AbsOmegaPiLhLS> GOmegaPiIndividualLS::omegaPiLhPtr() const {
      return _omegaPiLhPtr;
    }

    /********************************************************************************************/
    /**
     * A standard assignment operator
     *
     * @param cp A copy of another GOmegaPiIndividualLS object
     * @return A constant reference to this object
     */
    const GOmegaPiIndividualLS& GOmegaPiIndividualLS::operator=(const GOmegaPiIndividualLS& cp) {
      GOmegaPiIndividualLS::load_(&cp);
      return *this;
    }
  
    /********************************************************************************************/
    /**
     * Retrieval of fit parameters
     */
    bool GOmegaPiIndividualLS::getFitParams(OmegaPiDataLS::fitParamVal& fitParmVal) { 
      std::vector<double> theParms;
      // Extract the GDoubleCollection object. In a realistic scenario, you might want
      // to add error checks here upon first invocation.
      boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> vC = at<Gem::Geneva::GConstrainedDoubleObjectCollection>(0);
    
    
      std::vector<double> par;
      for(std::size_t i=0; i<vC->size(); i++){
	double value = vC->at(i)->value();
	par.push_back(value);
      }

      _omegaPiLhPtr->getFitParamVal(fitParmVal, par);
    
    
      return true;
    }

    /********************************************************************************************/
    /*
     * Printing of fit-parameters
     */
    void GOmegaPiIndividualLS::printFitParams(OmegaPiDataLS::fitParamVal& fitParmVal) {

      _omegaPiLhPtr->printFitParams(std::cout, fitParmVal);
      std::cout << std::endl;
      return;
    }
  
    /********************************************************************************************/
    /**
     * Loads the data of another GOmegaPiIndividualLS, camouflaged as a GObject.
     *
     * @param cp A copy of another GOmegaPiIndividualLS, camouflaged as a GObject
     */
    void GOmegaPiIndividualLS::load_(const Gem::Geneva::GObject* cp)
    {
      // Check that we are not accidently assigning this object to itself
      selfAssignmentCheck<GOmegaPiIndividualLS>(cp);
    
      // Load our parent's data
      GParameterSet::load_(cp);
    
      // Note: We do not need to load _omegaPiLhPtr here, as it has been initialized during
      // the construction of this object already.
    }
  
    /********************************************************************************************/
    /**
     * Creates a deep clone of this object
     *
     * @return A deep clone of this object, camouflaged as a GObject
     */
    Gem::Geneva::GObject* GOmegaPiIndividualLS::clone_() const {
      return new GOmegaPiIndividualLS(*this);
    }
  
    /********************************************************************************************/
    /**
     * The actual fitness calculation takes place here.
     *
     * @return The value of this object
     */
    double GOmegaPiIndividualLS::fitnessCalculation(){
      double result = 0.;
    
      OmegaPiDataLS::fitParamVal theFitParmValTmp;
      bool fitParamsSet=getFitParams(theFitParmValTmp);
      if (!fitParamsSet){
	std::cout << "fit parameters could not set properly " << std::endl;
	exit(0);
      }
      result=_omegaPiLhPtr->calcLogLh(theFitParmValTmp);
//       std::cout << "**** current fit params ***** " << std::endl;
//       printFitParams(theFitParmValTmp);
//       std::cout << std::endl;
      return result;
    }
  
    /********************************************************************************************/
    /**
     * The default constructor. Intentionally private and empty, as it is only needed for
     * serialization purposes.
     */
    GOmegaPiIndividualLS::GOmegaPiIndividualLS() :GParameterSet()
    {	/* nothing */ }

    /********************************************************************************************/
  
  } /* namespace Pawian */
} /* namespace Gem */
