#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/EtacToapi0Fit/EtacToapi0EventList.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
#include "ErrLogger/ErrLineLog.hh"


EtacToapi0Lh::EtacToapi0Lh(const EtacToapi0EventList* theEvtList) :
  _etacToapi0EventList(theEvtList)
{
}

EtacToapi0Lh::EtacToapi0Lh(const EtacToapi0Lh* theEtacToapi0Lh) :
  _etacToapi0EventList(theEtacToapi0Lh->getEventList())
{
}

EtacToapi0Lh::EtacToapi0Lh(boost::shared_ptr<EtacToapi0Lh> theEtacToapi0LhPtr):
  _etacToapi0EventList(theEtacToapi0LhPtr->getEventList())
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

  std::cout << "theParamVal.aMass= " << theParamVal.aMass << "\n"
	    << "theParamVal.aWidth= " << theParamVal.aWidth << "\n"
	    << "theParamVal.cont0spin= " << theParamVal.cont0spin << "\n"
	    << "theParamVal.cont1spin= " << theParamVal.cont1spin << "\n"
	    << "theParamVal.cont2spin= " << theParamVal.cont2spin << "\n";

  const std::vector<evt4Vec> data4Vecs=_etacToapi0EventList->getDataVecs();
  for (size_t it=0; it<data4Vecs.size(); it++)
    { 
      double intensity=calcEvtIntensity(data4Vecs[it], theParamVal);
      if (intensity>0.) logLH_data+=log10(intensity);
    }

  ErrMsg(debugging) << "logLH_data= " << logLH_data << endmsg;


  double LH_mc=0.;
  const std::vector<evt4Vec> mc4Vecs=_etacToapi0EventList->getMcVecs();
  for (size_t it=0; it<mc4Vecs.size(); it++)
    { 
      double intensity=calcEvtIntensity(mc4Vecs[it], theParamVal);
      LH_mc+=intensity;
    }
  ErrMsg(debugging) << "LH_mc= " << LH_mc << endmsg;

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log10(LH_mc/mc4Vecs.size());


  logLH=data4Vecs.size()/2.*(LH_mc/mc4Vecs.size()-1)*(LH_mc/mc4Vecs.size()-1)
    -logLH_data
    +data4Vecs.size()*logLH_mc_Norm;

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

