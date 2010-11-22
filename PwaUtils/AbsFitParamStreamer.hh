#ifndef AbsFitParamStreamer_H
#define AbsFitParamStreamer_H

#include <iostream>
#include <vector>
#include <map>
// #include <pair>
#include <string>
#include <boost/shared_ptr.hpp>
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

using namespace std;

typedef map<string, pair<double, double> > StringPairMap;

class AbsFitParamStreamer
{
  public:
  AbsFitParamStreamer(std::string& filePath);
  virtual ~AbsFitParamStreamer();

  protected:

  StringPairMap _stringPairMap;

  virtual void fillParamMap() = 0;
  virtual void StringParameterMap(const string &theName, const double firstVal, const double secondVal);
  virtual void fillParamMapAmps(std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSs, std::string& suffix, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& toFill);
  virtual void fillParamMapMass(std::string& name, pair<double, double>& toFill);
  virtual void fillParamFlatte(std::string& name, double& toFill);

  private:

};

#endif 
