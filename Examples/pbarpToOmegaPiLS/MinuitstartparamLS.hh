#ifndef MINUITSTARTPARAMLS_H
#define MINUITSTARTPARAMLS_H

#include "Examples/pbarpToOmegaPiLS/StreamparserLS.hh"
#include "Examples/pbarpToOmegaPiLS/FitparameterLS.hh"

class MinuitStartParamLS : public StreamParserLS, public FitParameterLS
{
  public:
    virtual void ParseStream(istream &theistream);
    
    inline void addMinuitUserParInMap(const string &theName, const double dMag, const double dPhi)
    {
      vector<double> tmpVec(2);
      tmpVec[0] = dMag;
      tmpVec[1] = dPhi;
      m_ParamMap[theName]=tmpVec;
    }

};

#endif // MINUITSTARTPARAM_H
