#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/OmegaPiLh.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"


OmegaPiLh::OmegaPiLh(boost::shared_ptr<const OmegaPiEventList> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0States> theStates) :
  _omegaPiEventListPtr(theEvtList),
  _omegaPi0StatesPtr(theStates)
{
}

OmegaPiLh::OmegaPiLh(boost::shared_ptr<OmegaPiLh> theOmegaPiLhPtr):
  _omegaPiEventListPtr(theOmegaPiLhPtr->getEventList()),
  _omegaPi0StatesPtr(theOmegaPiLhPtr->omegaPi0States())
{
}

OmegaPiLh::~OmegaPiLh()
{
}

double OmegaPiLh::calcLogLh(const OmegaPiData::fitParamVal& theParamVal){
  std::vector<OmegaPiData::OmPiEvtData> theData=_omegaPiEventListPtr->getDataVecs();
  std::vector<OmegaPiData::OmPiEvtData> theMCs=_omegaPiEventListPtr->getMcVecs();
 
  double logLH=0.;
  double logLH_data=0.;

  std::vector<OmegaPiData::OmPiEvtData>::iterator iterd;
  for (iterd=theData.begin(); iterd!=theData.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log10(intensity);
  } 

  double LH_mc=0.;
  
  std::vector<OmegaPiData::OmPiEvtData>::iterator iterm;
  for (iterm=theMCs.begin(); iterm!=theMCs.end(); ++iterm){
           double intensity=calcEvtIntensity((*iterm), theParamVal);
           LH_mc+=intensity;
         }

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log10(LH_mc/theMCs.size());

  logLH=theData.size()/2.*(LH_mc/theMCs.size()-1)*(LH_mc/theMCs.size()-1)
    -logLH_data
    +theData.size()*logLH_mc_Norm;

 return logLH;

}

double OmegaPiLh::calcEvtIntensity(OmegaPiData::OmPiEvtData& theData, const OmegaPiData::fitParamVal& theParamVal){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=theParamVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=theParamVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=theParamVal.omegaProdTriplet1;

  Spin lamgamma=-1;
  complex<double> singletAmpGM1=calcCoherentAmp(lamgamma,0, fitParmSinglet, theData);  
  complex<double> triplet0AmpGM1=calcCoherentAmp(lamgamma,0, fitParmTriplet0, theData);
  complex<double> tripletP1AmpGM1=calcCoherentAmp(lamgamma,1, fitParmTriplet1, theData);
  complex<double> tripletM1AmpGM1=calcCoherentAmp(lamgamma,-1, fitParmTriplet1, theData);

  lamgamma=1;
  complex<double> singletAmpGP1=calcCoherentAmp(lamgamma,0, fitParmSinglet, theData);  
  complex<double> triplet0AmpGP1=calcCoherentAmp(lamgamma,0, fitParmTriplet0, theData);
  complex<double> tripletP1AmpGP1=calcCoherentAmp(lamgamma,1, fitParmTriplet1, theData);
  complex<double> tripletM1AmpGP1=calcCoherentAmp(lamgamma,-1, fitParmTriplet1, theData);

  double result=2*norm(singletAmpGM1)+2*norm(triplet0AmpGM1)+norm(tripletP1AmpGM1)+norm(tripletM1AmpGM1);
  result+=2*norm(singletAmpGP1)+2*norm(triplet0AmpGP1)+norm(tripletP1AmpGP1)+norm(tripletM1AmpGP1);

  return result;  
}

complex<double> OmegaPiLh::calcCoherentAmp(Spin lamgamma, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData& theData){

  complex<double> result(0.,0.);
  for (Spin lamomega=-1; lamomega<=1; lamomega++){

    complex<double> omegaDecAmp=Clebsch(1,0,1,lamgamma,1, lamgamma)*Clebsch(1,lamgamma,0,0,1, lamgamma)*theData.Dfp[1][0][lamgamma];

    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
    for ( it=fitParm.begin(); it!=fitParm.end(); ++it){
      boost::shared_ptr<const JPCLS> theJPCLS=it->first;
      double theMag=it->second.first;
      double thePhi=it->second.second;
      complex<double> expiphi(cos(thePhi), sin(thePhi));
      complex<double> tmpResult=theMag*expiphi*Clebsch(theJPCLS->L,0,1, lamomega,theJPCLS->J, lamomega)*Clebsch(1,lamomega,0,0,1,lamomega)*theData.Dfp[theJPCLS->J][Minit][lamomega];     
      result+=tmpResult;                         
    } 
    result*=omegaDecAmp;   
  }

  return result;
}



void OmegaPiLh::print(std::ostream& os) const{
  os << "OmegaPiLh::print\n";
}
