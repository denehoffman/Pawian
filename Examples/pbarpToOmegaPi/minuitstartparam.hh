#ifndef MINUITSTARTPARAM_H
#define MINUITSTARTPARAM_H

#include "Examples/pbarpToOmegaPi/streamparser.hh"
#include "Examples/pbarpToOmegaPi/fitparameter.hh"

class minuitStartParam : public streamParser, public FitParameter
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
