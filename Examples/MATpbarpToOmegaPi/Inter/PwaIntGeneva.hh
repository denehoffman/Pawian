#ifndef _PWAINTGENEVA_H
#define _PWAINTGENEVA_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "Examples/MATpbarpToOmegaPi/GOmegaPiIndividual.hh"
#include "Examples/MATpbarpToOmegaPi/OmegaPiData.hh"

using namespace std;

class PwaIntGeneva {

public:
  /// Default Constructor (0x0)
  PwaIntGeneva();
  OmegaPiData::fitParamVal exec(vector<boost::shared_ptr<GOmegaPiIndividual> > in_fcn); 

  /** Destructor */
  virtual ~PwaIntGeneva();

 protected:

 private:
 // MOmegaPiFcn* fcn;
 // vector<string> paramNames;
};




#endif /* _PwaIntGeneva_H */
