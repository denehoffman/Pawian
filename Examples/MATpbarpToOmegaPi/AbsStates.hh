#ifndef _AbsStates_H
#define _AbsStates_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/text_iarchive.hpp> 
#include <boost/serialization/vector.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/MATpbarpToOmegaPi/serDataUtils.hh"


class AbsStates {

  friend class boost::serialization::access; 

  template <typename Archive> 
  void serialize(Archive & ar, const unsigned int) {
    using boost::serialization::make_nvp;

    /* 	 std::vector< boost::shared_ptr<const JPCLSM> > _allStates;
 	 std::vector< boost::shared_ptr<const JPCSM> > _allJPCSM;
 	 std::vector< boost::shared_ptr<const JPCLS> > _allJPCLS;
 	 std::vector< boost::shared_ptr<const jpcRes> >  _jpcStates;*/

    ar & make_nvp("allStates",_allStates);
    ar & make_nvp("allJPCSM",_allJPCSM);
    ar & make_nvp("allJPCLS",_allJPCLS);
    ar & make_nvp("jpcStates",_jpcStates);

    /* Add your own class-variables here in the following way:
      ar & BOOST_SERIALIZATION_NVP(myVar);
      or
      ar & make_nvp("myVar", myVar); // The latter form can be necessary when dealing with templates
    */
  }

public:

  /// Default Constructor
  AbsStates();
  AbsStates(int jmax);

  /** Destructor */
  virtual ~AbsStates();

  virtual void print(std::ostream& os) const;

  virtual std::vector< boost::shared_ptr<const JPCLSM> > allStates() const {return _allStates;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > jpcStates() const {return _jpcStates;}; 
  virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsStates() const {return _allJPCLS;}; 
  virtual std::vector< boost::shared_ptr<const JPCLSM> > extractJPCLSMStates(std::vector< boost::shared_ptr<const jpcRes> >&) const ;
  virtual std::vector< boost::shared_ptr<const JPCSM> > extractJPCSMStates(std::vector< boost::shared_ptr<const jpcRes> >&) const ;
  virtual std::vector< boost::shared_ptr<const JPCLS> > extractJPCLSStates(boost::shared_ptr<const jpcRes>) const;

 protected:
  std::vector< boost::shared_ptr<const JPCLSM> > _allStates;
  std::vector< boost::shared_ptr<const JPCSM> > _allJPCSM;
  std::vector< boost::shared_ptr<const JPCLS> > _allJPCLS;
  std::vector< boost::shared_ptr<const jpcRes> >  _jpcStates;

  virtual bool calcJPCs()=0;

 private:


};




#endif /* _AbsStates_H */
