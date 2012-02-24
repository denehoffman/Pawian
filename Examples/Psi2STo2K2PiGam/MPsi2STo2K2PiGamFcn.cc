//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#include "Minuit2/MnUserParameters.h"

#include "Examples/Psi2STo2K2PiGam/MPsi2STo2K2PiGamFcn.hh"
#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

MPsi2STo2K2PiGamFcn::MPsi2STo2K2PiGamFcn(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> psi2STo2K2PiGamLh) :
  _psi2STo2K2PiGamLhPtr(psi2STo2K2PiGamLh),
  _fcnCounter(0)
{
   if (0==_psi2STo2K2PiGamLhPtr) { Alert << "AbsPsi2STo2K2PiGamLh* _psi2STo2K2PiGamLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
//    (*_fcnCounter)=0;  
}

MPsi2STo2K2PiGamFcn::~MPsi2STo2K2PiGamFcn()
{
}

double MPsi2STo2K2PiGamFcn::operator()(const std::vector<double>& par) const
{

  _fcnCounter++;
  param2K2PiGam theFitParmValTmp;

  mutex1.lock(); 
  boost::shared_ptr<AbsPsi2STo2K2PiGamLh> tmpAbsLh=boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(_psi2STo2K2PiGamLhPtr->clone_());
  tmpAbsLh->setFitParamVal(theFitParmValTmp, par);    
  mutex1.unlock(); 

  double result=tmpAbsLh->calcLogLh(theFitParmValTmp);


  if (   _fcnCounter==1 || _fcnCounter%10 == 0) {
    mutex1.lock();
    DebugMsg << "logLh= " << result <<endmsg;  
    tmpAbsLh->printCurrentFitResult(theFitParmValTmp); 
    _psi2STo2K2PiGamLhPtr.reset();
    _psi2STo2K2PiGamLhPtr=tmpAbsLh;
    mutex1.unlock();  
  }

  if (  _fcnCounter%100 == 0) { 
   mutex1.lock(); 
    std::ofstream theStream ( "currentResult.dat");
    std::string theSuffix="Val"; 
    tmpAbsLh->dumpCurrentResult(theStream, theFitParmValTmp, theSuffix);
    mutex1.unlock();  
  }


//   _psi2STo2K2PiGamLhPtr.reset(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(_psi2STo2K2PiGamLhPtr->clone_()));
//   _psi2STo2K2PiGamLhPtr=boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(tmpAbsLh);
//   delete tmpAbsLh;

//   _psi2STo2K2PiGamLhPtr->setFitParamVal(theFitParmValTmp, par);
 
//   double result=_psi2STo2K2PiGamLhPtr->calcLogLh(theFitParmValTmp);

//   DebugMsg << "logLh= " << result <<endmsg;  
//   _psi2STo2K2PiGamLhPtr->printCurrentFitResult(theFitParmValTmp);

//   if (  (*_fcnCounter)%100 == 0) {
//     std::ofstream theStream ( "currentResult.dat");
//     std::string theSuffix="Val"; 
//     _psi2STo2K2PiGamLhPtr->dumpCurrentResult(theStream, theFitParmValTmp, theSuffix);
//   }
  return result;
}

double MPsi2STo2K2PiGamFcn::Up() const 
{
return .5;
}

