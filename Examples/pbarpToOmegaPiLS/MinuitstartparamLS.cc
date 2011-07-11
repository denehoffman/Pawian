#include "Examples/pbarpToOmegaPiLS/MinuitstartparamLS.hh"
#include "ErrLogger/ErrLogger.hh"

void MinuitStartParamLS::ParseStream(istream &theistream)
{
      string strTmp, strState;
      bool bError=false;
      double dMag=0.0, dPhi=0.0;
  do
  {
    theistream >> strState;
    if (theistream.eof()) break;
    if (strState.empty()) bError=true;
    theistream >> strTmp;
    dMag = atof(strTmp.c_str());
    if (theistream.eof()) break;
    theistream >> strTmp;
    dPhi = atof(strTmp.c_str());
    if (theistream.eof()) break;
    if (!bError) 
    {
      addMinuitUserParInMap(strState,dMag,dPhi);
    }
  }
  while(!theistream.eof());
}

