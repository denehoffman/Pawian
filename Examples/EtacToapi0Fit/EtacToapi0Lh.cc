#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/EtacToapi0Fit/EtacToapi0EventList.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
#include "ErrLogger/ErrLineLog.hh"


EtacToapi0Lh::EtacToapi0Lh(boost::shared_ptr<const EtacToapi0EventList> theEvtList) :
  _etacToapi0EventListPtr(theEvtList)
{
}

EtacToapi0Lh::EtacToapi0Lh(boost::shared_ptr<EtacToapi0Lh> theEtacToapi0LhPtr):
  _etacToapi0EventListPtr(theEtacToapi0LhPtr->getEventList())
{
}

EtacToapi0Lh::~EtacToapi0Lh()
{
}

double EtacToapi0Lh::calcEvtIntensity(const evt4Vec& theEvtVec, const fitParamVal& theParamVal)
{
  double spin2cont=theParamVal.cont2spin;
  double spin1cont=theParamVal.cont1spin;
  double spin0cont=theParamVal.cont0spin;
  
  double intensity=0;

  Tensor<complex<double> > amp(0); // rank-0 Tensor for the amplitude

  amp = theEvtVec.spin2_0_amp*BreitWigner(theEvtVec.a2_0_4Vec, theParamVal.aMass, theParamVal.aWidth);
  intensity += spin2cont*norm(amp(0,0));

  amp = theEvtVec.spin2_1_amp*BreitWigner(theEvtVec.a2_1_4Vec, theParamVal.aMass, theParamVal.aWidth);
  intensity += spin2cont*norm(amp(0,0));

  amp = theEvtVec.spin1_0_amp*BreitWigner(theEvtVec.a2_0_4Vec, theParamVal.aMass, theParamVal.aWidth);
  intensity += spin1cont*norm(amp(0,0));

  amp = theEvtVec.spin1_1_amp*BreitWigner(theEvtVec.a2_1_4Vec, theParamVal.aMass, theParamVal.aWidth);
  intensity += spin1cont*norm(amp(0,0));

  //dont know how to normalize
  amp = 0.1*BreitWigner(theEvtVec.a2_0_4Vec, theParamVal.aMass, theParamVal.aWidth);
  intensity += spin0cont*norm(amp(0,0));

  amp = 0.1*BreitWigner(theEvtVec.a2_1_4Vec, theParamVal.aMass, theParamVal.aWidth);
  intensity += spin0cont*norm(amp(0,0));


  return intensity;
}



double EtacToapi0Lh::calcLogLh(const fitParamVal& theParamVal)
{
  double logLH=0.;
  double logLH_data=0.;


  const std::vector<evt4Vec> data4Vecs=_etacToapi0EventListPtr->getDataVecs();


  std::vector<evt4Vec>::const_iterator iterd;
  for (iterd=data4Vecs.begin(); iterd!=data4Vecs.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log10(intensity);
  }  

  double LH_mc=0.;
  const std::vector<evt4Vec> mc4Vecs=_etacToapi0EventListPtr->getMcVecs();
  
  std::vector<evt4Vec>::const_iterator iterm;
  for (iterm=mc4Vecs.begin(); iterm!=mc4Vecs.end(); ++iterm){
	   double intensity=calcEvtIntensity((*iterm), theParamVal);
	   LH_mc+=intensity;
	 }


  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log10(LH_mc/mc4Vecs.size());


  logLH=data4Vecs.size()/2.*(LH_mc/mc4Vecs.size()-1)*(LH_mc/mc4Vecs.size()-1)
    -logLH_data
    +data4Vecs.size()*logLH_mc_Norm;

  ErrMsg(debugging) << "theParamVal.aMass= " << theParamVal.aMass << "\n"
		    << "theParamVal.aWidth= " << theParamVal.aWidth << "\n"
		    << "theParamVal.cont0spin= " << theParamVal.cont0spin << "\n"
		    << "theParamVal.cont1spin= " << theParamVal.cont1spin << "\n"
		    << "theParamVal.cont2spin= " << theParamVal.cont2spin << "\n"
		    << "logLH_data= " << logLH_data << "\n"  
		    << "LH_mc= " << LH_mc << "\n"
		    << "lgLh= " << logLH << endmsg;
 return logLH;
}


bool EtacToapi0Lh::setFitParamVal(fitParamVal& fitParamVal, const std::vector<double>& par)
{
  bool result=true;
  fitParamVal.aMass=par[0];
  fitParamVal.aWidth=par[1];
  fitParamVal.cont0spin=par[2];
  fitParamVal.cont1spin=par[3];
  fitParamVal.cont2spin=par[4];
  return result;
} 

