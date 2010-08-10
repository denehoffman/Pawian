#include <string>
#include <vector>
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "Examples/MATpbarpToOmegaPi/MOmegaPiFcn.hh"
#include "Examples/MATpbarpToOmegaPi/PwaIntMinuit.hh"

using namespace ROOT::Minuit2;
using namespace std;

//TODO: 
//general fcn class

PwaIntMinuit::PwaIntMinuit(){
 // for(int i=in_paramNames.size()-1; i>-1; i--)
  //  paramNames.push_back(in_paramNames.at(i));
}

PwaIntMinuit::~PwaIntMinuit()
{
}

const vector<double> PwaIntMinuit::exec(MOmegaPiFcn* in_fcn){
  
  MnUserParameters upar;
  in_fcn->setMnUsrParams(upar);

  MnMigrad migrad(*in_fcn, upar);
 // Info <<"start migrad "<< endmsg;
  FunctionMinimum min = migrad();

 if(!min.IsValid()) {
   //try with higher strategy
 //  Info <<"FM is invalid, try with strategy = 2."<< endmsg;
   MnMigrad migrad2(*in_fcn, min.UserState(), MnStrategy(2));
   min = migrad2();
 }

 MnUserParameters finalUsrParameters=min.UserParameters();
 
 return finalUsrParameters.Params();
}
