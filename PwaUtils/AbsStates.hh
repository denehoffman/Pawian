#ifndef _AbsStates_H
#define _AbsStates_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"


class AbsStates {

public:

  /// Default Constructor
  AbsStates();
  AbsStates(int jmax);

  /** Destructor */
  virtual ~AbsStates();

  virtual void print(std::ostream& os) const;

  virtual std::vector< boost::shared_ptr<const JPCSML> > allStates() {return _allStates;};
  virtual std::vector< boost::shared_ptr<const JPCSML> > singletStates() {return _singletStates;};
  virtual std::vector< boost::shared_ptr<const JPCSML> > tripletM0States() {return _tripletM0States;};
  virtual std::vector< boost::shared_ptr<const JPCSML> > tripletMp1States() {return _tripletMp1States;};
  virtual std::vector< boost::shared_ptr<const JPCSML> > tripletMm1States() {return _tripletMm1States;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > jpcStates() {return _jpcStates;}; 
  virtual std::vector< boost::shared_ptr<const JPCSML> > extractJPCSMLStates(std::vector< boost::shared_ptr<const jpcRes> >&);
  virtual std::vector< boost::shared_ptr<const JPCSM> > extractJPCSMStates(std::vector< boost::shared_ptr<const jpcRes> >&);

 protected:
  std::vector< boost::shared_ptr<const JPCSML> > _allStates;
  std::vector< boost::shared_ptr<const JPCSML> > _singletStates;
  std::vector< boost::shared_ptr<const JPCSML> > _tripletM0States;
  std::vector< boost::shared_ptr<const JPCSML> > _tripletMp1States;
  std::vector< boost::shared_ptr<const JPCSML> > _tripletMm1States;
  std::vector< boost::shared_ptr<const JPCSM> > _allJPCSM;
  std::vector< boost::shared_ptr<const jpcRes> >  _jpcStates;

  virtual bool calcJPCs()=0;

 private:


};




#endif /* _AbsStates_H */
