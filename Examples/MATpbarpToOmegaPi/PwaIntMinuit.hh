#ifndef _PWAINTMINUIT_H
#define _PWAINTMINUIT_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "Examples/MATpbarpToOmegaPi/MOmegaPiFcn.hh"

using namespace ROOT::Minuit2;
using namespace std;

class PwaIntMinuit {

public:
  /// Default Constructor (0x0)
  PwaIntMinuit();
  const vector<double> exec(MOmegaPiFcn* in_fcn); 

  /** Destructor */
  virtual ~PwaIntMinuit();

 protected:

 private:
 // MOmegaPiFcn* fcn;
 // vector<string> paramNames;
};




#endif /* _PwaIntMinuit_H */
