#include <getopt.h>

#include "Examples/pbarpToOmegaPiLS/OmegaTo3PiLhPi0GammaLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "Examples/pbarpToOmegaPiLS/pbarpToOmegaPi0StatesLS.hh"
#include "ErrLogger/ErrLogger.hh"

// #include <geneva/GConstrainedDoubleObject.hpp>


OmegaTo3PiLhPi0GammaLS::OmegaTo3PiLhPi0GammaLS(boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0StatesLS> theStates) :
  AbsOmegaPiLhLS(theEvtList, theStates)
{
}

OmegaTo3PiLhPi0GammaLS::OmegaTo3PiLhPi0GammaLS(boost::shared_ptr<OmegaTo3PiLhPi0GammaLS> theOmegaTo3PiLhPi0GammaLSPtr):
  AbsOmegaPiLhLS(theOmegaTo3PiLhPi0GammaLSPtr)
{
}

OmegaTo3PiLhPi0GammaLS::~OmegaTo3PiLhPi0GammaLS()
{
}

double OmegaTo3PiLhPi0GammaLS::calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal){

  double result=AbsOmegaPiLhLS::calcLogLh(theParamVal);
  if (_globalItCounter%10 == 0) printFitParams(std::cout, theParamVal);

//   if (_globalItCounter%1000 == 0){
//     std::ofstream theStream ("currentResult.dat");
//     dumpCurrentResult(theStream, theParamVal);
//     theStream.close();
//   }

  return result; 
}


double OmegaTo3PiLhPi0GammaLS::calcEvtIntensity(OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal){

  Spin lamOmegaDecParticle=0;
  complex<double> singletAmp=calcCoherentAmp(lamOmegaDecParticle,0, theParamVal, _singlet_JPCLS_States, theData);
  complex<double> triplet0Amp=calcCoherentAmp(lamOmegaDecParticle,0, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> tripletP1Amp=calcCoherentAmp(lamOmegaDecParticle,1, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> tripletM1Amp=calcCoherentAmp(lamOmegaDecParticle,-1, theParamVal, _tripletm1_JPCLS_States, theData);


  double result=norm(singletAmp)+norm(triplet0Amp)+norm(tripletP1Amp)+norm(tripletM1Amp);  
  return result;  
}


complex<double> OmegaTo3PiLhPi0GammaLS::calcCoherentAmp(Spin lamOmDec, Spin Minit, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData){

  complex<double> result(0.,0.);

  for (Spin lamomega=-1; lamomega<=1; lamomega++){

    complex<double> omegaDecAmp=conj(theData->Dfd[1][lamomega][0]);// Clebsch(1, 0,0,0,1,0)=1;Clebsch(0, 0,  0, 0, 0, 0)=1

    std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator it;
     for ( it=theJPCLSlsStates.begin(); it!=theJPCLSlsStates.end(); ++it){


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


void OmegaTo3PiLhPi0GammaLS::print(std::ostream& os) const{
  os << "OmegaTo3PiLhPi0GammaLS\n";
}
