#ifndef _AbsLh_H
#define _AbsLh_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"


class AbsLh {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsLh(boost::shared_ptr<const EvtDataBaseList>);
  AbsLh(boost::shared_ptr<AbsLh>);

  /** Destructor */
  virtual ~AbsLh();

  virtual AbsLh* clone_() const =0;


  // Getters:
  
  virtual double calcLogLh(fitParams& theParamVal);

  virtual double calcEvtIntensity(EvtData* theData, fitParams& theParamVal)=0;

  virtual boost::shared_ptr<const EvtDataBaseList> getEventList() const {return _evtListPtr;}

  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr)=0;
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
