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

  /** Destructor */
  virtual ~AbsStates();
  virtual void print(std::ostream& os) const;

//   virtual std::vector< boost::shared_ptr<const JPCLSM> > allStates() const {return _allStates;};
//   virtual std::vector< boost::shared_ptr<const jpcRes> > jpcStates() const {return _jpcStates;}; 
//   virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsStates() const {return _allJPCLS;}; 
//   virtual std::vector< boost::shared_ptr<const JPCLSM> > extractJPCLSMStates(std::vector< boost::shared_ptr<const jpcRes> >&) const ;
//    virtual std::vector< boost::shared_ptr<const JPCSM> > extractJPCSMStates(std::vector< boost::shared_ptr<const jpcRes> >&) const ;
//   virtual std::vector< boost::shared_ptr<const JPCLS> > extractJPCLSStates(boost::shared_ptr<const jpcRes>) const;

 protected:
  virtual bool calcStates()=0;
//   std::vector< boost::shared_ptr<const JPCLSM> > _allStates;
//   std::vector< boost::shared_ptr<const JPCSM> > _allJPCSM;
//   std::vector< boost::shared_ptr<const JPCLS> > _allJPCLS;
//   std::vector< boost::shared_ptr<const jpcRes> >  _jpcStates;


 private:


};




#endif /* _AbsStates_H */
