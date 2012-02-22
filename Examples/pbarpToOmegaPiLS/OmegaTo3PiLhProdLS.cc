#include <getopt.h>

#include "Examples/pbarpToOmegaPiLS/OmegaTo3PiLhProdLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "Examples/pbarpToOmegaPiLS/pbarpToOmegaPi0StatesLS.hh"
#include "ErrLogger/ErrLogger.hh"

// #include <geneva/GConstrainedDoubleObject.hpp>


OmegaTo3PiLhProdLS::OmegaTo3PiLhProdLS(boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0StatesLS> theStates) :
  AbsOmegaPiLhLS(theEvtList, theStates, "OmegaTo3PiLhProd")
{
}

OmegaTo3PiLhProdLS::OmegaTo3PiLhProdLS(boost::shared_ptr<OmegaTo3PiLhProdLS> theOmegaTo3PiLhProdLSPtr):
  AbsOmegaPiLhLS(theOmegaTo3PiLhProdLSPtr, "OmegaTo3PiLhProd")
{
}

OmegaTo3PiLhProdLS::~OmegaTo3PiLhProdLS()
{
}

double OmegaTo3PiLhProdLS::calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal){

  double result=AbsOmegaPiLhLS::calcLogLh(theParamVal);
  if (_globalItCounter%10 == 0) printFitParams(std::cout, theParamVal);

//   if (_globalItCounter%1000 == 0){
//     std::ofstream theStream ("currentResult.dat");
//     dumpCurrentResult(theStream, theParamVal);
//     theStream.close();
//   }

  return result; 
}


double OmegaTo3PiLhProdLS::calcEvtIntensity(OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal){

  double result=0.;

  for (Spin lamOmega=-1; lamOmega<=1; ++lamOmega){
    result+= norm(calcCoherentAmp(lamOmega,0, theParamVal, _singlet_JPCLS_States, theData)); //singlet amp
    result+= norm(calcCoherentAmp(lamOmega,0, theParamVal, _triplet0_JPCLS_States, theData)); //triplet0
    result+= norm(calcCoherentAmp(lamOmega,1, theParamVal, _tripletp1_JPCLS_States, theData)); //triplet+1
    result+= norm(calcCoherentAmp(lamOmega,-1, theParamVal, _tripletm1_JPCLS_States, theData)); //triplet-1
  }

  return result;  
}


complex<double> OmegaTo3PiLhProdLS::calcCoherentAmp(Spin lamOm, Spin Minit, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData){

  complex<double> result(0.,0.);

  std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator it;
  for ( it=theJPCLSlsStates.begin(); it!=theJPCLSlsStates.end(); ++it){
    
    
    std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator itmap;
    itmap = theParamVal.lsParam.find((*it));
    pair<double, double> fitPair= (*itmap).second; 
    
    if (fabs(lamOm)> (*it)->J ) continue;
    
    complex<double> omegaProdAmp=calcOmegaProdPartAmp(Minit, lamOm, (*it), fitPair, theData);
    
    result += omegaProdAmp;
    
  }

  

  return result;
  }


void OmegaTo3PiLhProdLS::print(std::ostream& os) const{
  os << "OmegaTo3PiLhProdLS\n";
}
