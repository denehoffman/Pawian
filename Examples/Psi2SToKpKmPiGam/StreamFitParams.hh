#ifndef StreamFitParams_H
#define StreamFitParams_H

#include <iostream>
#include <vector>
#include <map>
//#include <pair>
#include <string>
#include <boost/shared_ptr.hpp>
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamStates.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "PwaUtils/DataUtils.hh"

using namespace std;

typedef map<string, pair<double, double> > StringPairMap;

class StreamFitParams
{
  public:
  StreamFitParams(std::string& filePath, boost::shared_ptr<const Psi2SToKpKmPiGamStates>);

  Psi2SToKpKmPiGamData::fitParamVal getFitParamVal() { return _paramVal;}
  Psi2SToKpKmPiGamData::fitParamVal getFitParamErr() { return _paramErr;}


  protected:

  private:
  Psi2SToKpKmPiGamData::fitParamVal _paramVal;
  Psi2SToKpKmPiGamData::fitParamVal _paramErr;
  StringPairMap _stringPairMap;

  void StringParameterMap(const string &theName, const double firstVal, const double secondVal);
  void fillParamMap(boost::shared_ptr<const Psi2SToKpKmPiGamStates>);
  void fillParamMapAmps(std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSs, std::string& suffix, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& toFill);
  void fillParamMapMass(std::string& name, pair<double, double>& toFill);
  void fillParamFlatte(std::string& name, double& toFill);
};

#endif // MINUITSTARTPARAM_H
