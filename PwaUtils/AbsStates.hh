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

  virtual std::vector< boost::shared_ptr<const JPCLSM> > allStates() {return _allStates;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > jpcStates() {return _jpcStates;}; 
  virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsStates() {return _allJPCLS;}; 
  virtual std::vector< boost::shared_ptr<const JPCLSM> > extractJPCLSMStates(std::vector< boost::shared_ptr<const jpcRes> >&);
  virtual std::vector< boost::shared_ptr<const JPCSM> > extractJPCSMStates(std::vector< boost::shared_ptr<const jpcRes> >&);
  virtual std::vector< boost::shared_ptr<const JPCLS> > extractJPCLSStates(boost::shared_ptr<const jpcRes>);

 protected:
  std::vector< boost::shared_ptr<const JPCLSM> > _allStates;
  std::vector< boost::shared_ptr<const JPCSM> > _allJPCSM;
  std::vector< boost::shared_ptr<const JPCLS> > _allJPCLS;
  std::vector< boost::shared_ptr<const jpcRes> >  _jpcStates;

  virtual bool calcJPCs()=0;

 private:


};




#endif /* _AbsStates_H */
