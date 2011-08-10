#ifndef _MIMINUIT_H
#define _MIMINUIT_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIData.hh"
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIBase.hh"
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIMinuitFcn.hh"

using namespace ROOT::Minuit2;
using namespace std;

class MIMinuit : public MIBase {

public:
  /// Default Constructor (0x0)
  MIMinuit(boost::shared_ptr<MIData> theData);  
  virtual const double exec(unsigned int Npar, double* par,  double* min, double* max, double* err); 

  /** Destructor */
  virtual ~MIMinuit();

 protected:

 private:
  MIMinuitFcn _myFcn;
 // vector<string> paramNames;
};

#endif /* _MIMinuit_H */
