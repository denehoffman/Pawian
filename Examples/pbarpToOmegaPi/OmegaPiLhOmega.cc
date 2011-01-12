#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/OmegaPiLhOmega.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "ErrLogger/ErrLogger.hh"

OmegaPiLhOmega::OmegaPiLhOmega(boost::shared_ptr<const OmegaPiEventList> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0States> theStates) :
  AbsOmegaPiLh(theEvtList, theStates)
{
}

OmegaPiLhOmega::OmegaPiLhOmega(boost::shared_ptr<OmegaPiLhOmega> theOmegaPiLhOmegaPtr):
  AbsOmegaPiLh(theOmegaPiLhOmegaPtr)
{
}

OmegaPiLhOmega::~OmegaPiLhOmega()
{
}


double OmegaPiLhOmega::calcEvtIntensity(OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=theParamVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=theParamVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=theParamVal.omegaProdTriplet1;


  Spin lamOmega=-1;
  complex<double> singletAmpOM1=calcCoherentAmp(lamOmega,0, fitParmSinglet, theData);
  complex<double> triplet0AmpOM1=calcCoherentAmp(lamOmega,0, fitParmTriplet0, theData);
  complex<double> tripletP1AmpOM1=calcCoherentAmp(lamOmega,1, fitParmTriplet1, theData);
  complex<double> tripletM1AmpOM1=calcCoherentAmp(lamOmega,-1, fitParmTriplet1, theData);

  lamOmega=0;
  complex<double> singletAmpO0=calcCoherentAmp(lamOmega,0, fitParmSinglet, theData);  
  complex<double> triplet0AmpO0=calcCoherentAmp(lamOmega,0, fitParmTriplet0, theData);
  complex<double> tripletP1AmpO0=calcCoherentAmp(lamOmega,1, fitParmTriplet1, theData);
  complex<double> tripletM1AmpO0=calcCoherentAmp(lamOmega,-1, fitParmTriplet1, theData);

  lamOmega=1;
  complex<double> singletAmpOP1=calcCoherentAmp(lamOmega,0, fitParmSinglet, theData);  
  complex<double> triplet0AmpOP1=calcCoherentAmp(lamOmega,0, fitParmTriplet0, theData);
  complex<double> tripletP1AmpOP1=calcCoherentAmp(lamOmega,1, fitParmTriplet1, theData);
  complex<double> tripletM1AmpOP1=calcCoherentAmp(lamOmega,-1, fitParmTriplet1, theData);


  double result=norm(singletAmpOM1)+norm(triplet0AmpOM1)+norm(tripletP1AmpOM1)+norm(tripletM1AmpOM1);
  result+=norm(singletAmpO0)+norm(triplet0AmpO0)+norm(tripletP1AmpO0)+norm(tripletM1AmpO0);  
  result+=norm(singletAmpOP1)+norm(triplet0AmpOP1)+norm(tripletP1AmpOP1)+norm(tripletM1AmpOP1);

  return result;  
}


complex<double> OmegaPiLhOmega::calcCoherentAmp(Spin lamomega, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData){

  complex<double> result(0.,0.);


    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
    
    for ( it=fitParm.begin(); it!=fitParm.end(); ++it){
      boost::shared_ptr<const JPCLS> theJPCLS=it->first;
      if (fabs(lamomega)>theJPCLS->J) continue;
      double theMag=it->second.first;
      double thePhi=it->second.second;
      complex<double> expiphi(cos(thePhi), sin(thePhi));
      result+=sqrt(2*theJPCLS->L+1)*theMag*expiphi*Clebsch(theJPCLS->L,0,1, lamomega,theJPCLS->J, lamomega)*conj(theData->Dfp[theJPCLS->J][Minit][lamomega]); //Clebsch(1,lamomega,0,0,1,lamomega)=1    
    }

  return result;
}


void OmegaPiLhOmega::print(std::ostream& os) const{
  os << "OmegaPiLhOmega::print\n";
}
