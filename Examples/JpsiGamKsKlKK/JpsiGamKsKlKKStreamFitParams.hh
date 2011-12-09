#ifndef JpsiGamKsKlKKStreamFitParams_H
#define JpsiGamKsKlKKStreamFitParams_H

#include <iostream>
#include <vector>
#include <map>
//#include <pair>
#include <string>
#include <boost/shared_ptr.hpp>
#include "PwaUtils/AbsFitParamStreamer.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKFitParams.hh"


class JpsiGamKsKlKKStreamFitParams : public AbsFitParamStreamer{
  
public:
  JpsiGamKsKlKKStreamFitParams(std::string&);
  virtual ~JpsiGamKsKlKKStreamFitParams();
  
  paramGamKsKlKK getFitParamVal() { return _paramVal;}
  paramGamKsKlKK getFitParamErr() { return _paramErr;}
  
protected:
  virtual void fillParamMap();

private:
  paramGamKsKlKK _paramVal;
  paramGamKsKlKK _paramErr;

  JpsiGamKsKlKKFitParams _fitparamsGamKsKlKK;
};

#endif 
