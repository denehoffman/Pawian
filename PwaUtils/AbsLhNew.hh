#ifndef _AbsLhNew_H
#define _AbsLhNew_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBaseNew.hh"


class AbsLhNew {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsLhNew(boost::shared_ptr<const EvtDataBaseList>);
  AbsLhNew(boost::shared_ptr<AbsLhNew>);

  /** Destructor */
  virtual ~AbsLhNew();

  virtual AbsLhNew* clone_() const =0;


  // Getters:
  
  virtual double calcLogLh(fitParamsNew& theParamVal);

  virtual double calcEvtIntensity(EvtData* theData, fitParamsNew& theParamVal)=0;

  virtual boost::shared_ptr<const EvtDataBaseList> getEventList() const {return _evtListPtr;}

  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr)=0;
  virtual void print(std::ostream& os) const=0;

protected:

  boost::shared_ptr<const EvtDataBaseList> _evtListPtr;

 
  std::vector<EvtData*> _evtDataVec;
  std::vector<EvtData*> _evtMCVec;
  std::map<const std::string, bool> _hypMap;

  virtual void setHyps( const std::map<const std::string, bool>& theMap, bool& theHyp, std::string& theKey);


private:


};

#endif

//  LocalWords:  ifndef
