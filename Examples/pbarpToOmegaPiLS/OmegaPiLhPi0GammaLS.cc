#include <getopt.h>

#include "Examples/pbarpToOmegaPiLS/OmegaPiLhPi0GammaLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "Examples/pbarpToOmegaPiLS/pbarpToOmegaPi0StatesLS.hh"
#include "ErrLogger/ErrLogger.hh"

// #include <geneva/GConstrainedDoubleObject.hpp>


OmegaPiLhPi0GammaLS::OmegaPiLhPi0GammaLS(boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0StatesLS> theStates) :
  AbsOmegaPiLhLS(theEvtList, theStates, "OmegaPiLhGamma")
{
}

OmegaPiLhPi0GammaLS::OmegaPiLhPi0GammaLS(boost::shared_ptr<OmegaPiLhPi0GammaLS> theOmegaPiLhPi0GammaLSPtr):
  AbsOmegaPiLhLS(theOmegaPiLhPi0GammaLSPtr, "OmegaPiLhGamma")
{
}

OmegaPiLhPi0GammaLS::~OmegaPiLhPi0GammaLS()
{
}

double OmegaPiLhPi0GammaLS::calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal){

  double result=AbsOmegaPiLhLS::calcLogLh(theParamVal);
  if (_globalItCounter%10 == 0) printFitParams(std::cout, theParamVal);

//   if (_globalItCounter%1000 == 0){
//     std::ofstream theStream ("currentResult.dat");
//     dumpCurrentResult(theStream, theParamVal);
//     theStream.close();
//   }

  return result; 
}


double OmegaPiLhPi0GammaLS::calcEvtIntensity(OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal){

  Spin lamgamma=-1;
  complex<double> singletAmpGM1=calcCoherentAmp(lamgamma,0, theParamVal, _singlet_JPCLS_States, theData);
  complex<double> triplet0AmpGM1=calcCoherentAmp(lamgamma,0, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> tripletP1AmpGM1=calcCoherentAmp(lamgamma,1, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> tripletM1AmpGM1=calcCoherentAmp(lamgamma,-1, theParamVal, _tripletm1_JPCLS_States, theData);

  lamgamma=1;
  complex<double> singletAmpGP1=calcCoherentAmp(lamgamma,0, theParamVal, _singlet_JPCLS_States, theData);  
  complex<double> triplet0AmpGP1=calcCoherentAmp(lamgamma,0, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> tripletP1AmpGP1=calcCoherentAmp(lamgamma,1, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> tripletM1AmpGP1=calcCoherentAmp(lamgamma,-1, theParamVal, _tripletm1_JPCLS_States, theData);

  double result=norm(singletAmpGM1)+norm(triplet0AmpGM1)+norm(tripletP1AmpGM1)+norm(tripletM1AmpGM1);
  
  result+=norm(singletAmpGP1)+norm(triplet0AmpGP1)+norm(tripletP1AmpGP1)+norm(tripletM1AmpGP1);

  return result;  
}


complex<double> OmegaPiLhPi0GammaLS::calcCoherentAmp(Spin lamgamma, Spin Minit, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData){

  complex<double> result(0.,0.);

  for (Spin lamomega=-1; lamomega<=1; lamomega++){

    complex<double> omegaDecAmp=Clebsch(1,0,1,lamgamma,1, lamgamma)*conj(theData->Dfd[1][lamomega][lamgamma]);// Clebsch(1,lamgamma,0,0,1, lamgamma)=1
    

    std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator it;
     for ( it=theJPCLSlsStates.begin(); it!=theJPCLSlsStates.end(); ++it){

//        pair<double, double> fitPair= theParamVal.lsParam[(*it)];

       std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator itmap;
       itmap = theParamVal.lsParam.find((*it));
       pair<double, double> fitPair= (*itmap).second; 
      
       if (fabs(lamomega)> (*it)->J ) continue;

       complex<double> omegaProdAmp=calcOmegaProdPartAmp(Minit, lamomega, (*it), fitPair, theData);
      
       result += omegaProdAmp*omegaDecAmp;
     }

  }

  return result;
  }


void OmegaPiLhPi0GammaLS::print(std::ostream& os) const{
  os << "OmegaPiLhPi0GammaLS\n";
}
