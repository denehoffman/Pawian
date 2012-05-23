#ifndef _AbsXdecAmp_H
#define _AbsXdecAmp_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBaseNew.hh"


class AbsXdecAmp {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsXdecAmp(const std::string& name, const std::vector<std::string>& hypVec);

  /** Destructor */
  virtual ~AbsXdecAmp();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, fitParamsNew& theParamVal)=0;
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr)=0;
  virtual void print(std::ostream& os) const=0;
  const std::string name() const {return _name;}

protected:

  const std::string _name;
  const std::vector<std::string> _hypVec;

private:


};

#endif


