#ifndef _AbsLhNew_H
#define _AbsLhNew_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"


class AbsLhNew {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsLhNew(boost::shared_ptr<const EvtDataBaseListNew>);
  AbsLhNew(boost::shared_ptr<AbsLhNew>);

  /** Destructor */
  virtual ~AbsLhNew();

  virtual AbsLhNew* clone_() const =0;


  // Getters:
  
  virtual double calcLogLh(fitParamsNew& theParamVal);

  virtual double calcEvtIntensity(EvtDataNew* theData, fitParamsNew& theParamVal)=0;

  virtual boost::shared_ptr<const EvtDataBaseListNew> getEventList() const {return _evtListPtr;}

  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr)=0;
  virtual void print(std::ostream& os) const=0;

protected:

  boost::shared_ptr<const EvtDataBaseListNew> _evtListPtr;

 
  std::vector<EvtDataNew*> _evtDataVec;
  std::vector<EvtDataNew*> _evtMCVec;
  std::map<const std::string, bool> _hypMap;

  virtual void setHyps( const std::map<const std::string, bool>& theMap, bool& theHyp, std::string& theKey);


private:


};

#endif

//  LocalWords:  ifndef
