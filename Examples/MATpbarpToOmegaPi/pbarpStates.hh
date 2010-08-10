#ifndef _pbarpStates_H
#define _pbarpStates_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/text_iarchive.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/MATpbarpToOmegaPi/AbsStates.hh"
#include "Examples/MATpbarpToOmegaPi/serDataUtils.hh"

#include "Examples/MATpbarpToOmegaPi/serSpin.hh"


class pbarpStates : public AbsStates {

  friend class boost::serialization::access;
  
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) {
    using boost::serialization::make_nvp;

    ar & BOOST_SERIALIZATION_NVP(_jmax);
    ar & BOOST_SERIALIZATION_NVP(_pbarJPC);
    ar & BOOST_SERIALIZATION_NVP(_pJPC);
    ar & BOOST_SERIALIZATION_NVP(_singletStates);
    ar & BOOST_SERIALIZATION_NVP(_tripletM0States);
    ar & BOOST_SERIALIZATION_NVP(_tripletMp1States);
    ar & BOOST_SERIALIZATION_NVP(_tripletMm1States);

    /* Add your own class-variables here in the following way:
      ar & BOOST_SERIALIZATION_NVP(myVar);
      or
      ar & make_nvp("myVar", myVar); // The latter form can be necessary when dealing with templates
    */
  }

public:

  /// Default Constructor
  pbarpStates();
  pbarpStates(int jmax);

  /** Destructor */
  virtual ~pbarpStates();

  virtual void print(std::ostream& os) const;

  virtual std::vector< boost::shared_ptr<const jpcRes> > singletStates() const {return _singletStates;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > tripletM0States() const {return _tripletM0States;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > tripletMp1States() const {return _tripletMp1States;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > tripletMm1States() const {return _tripletMm1States;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > extractSingletStates(std::vector< boost::shared_ptr<const jpcRes> >&) const; 
  virtual std::vector< boost::shared_ptr<const jpcRes> > extractTripletM0States(std::vector< boost::shared_ptr<const jpcRes> >&) const;
  virtual std::vector< boost::shared_ptr<const jpcRes> > extractTripletM1States(std::vector< boost::shared_ptr<const jpcRes> >&) const;

 protected:
  virtual bool calcJPCs();

 private:
  int _jmax;
  jpcRes _pbarJPC;
  jpcRes _pJPC;

  vector<LS> myGetValidLS(const serSpin &__j,int __parity,const serSpin &__s1,int __p1,
		      const serSpin &__s2,int __p2);

  std::vector< boost::shared_ptr<const jpcRes> > _singletStates;
  std::vector< boost::shared_ptr<const jpcRes> > _tripletM0States;
  std::vector< boost::shared_ptr<const jpcRes> > _tripletMp1States;
  std::vector< boost::shared_ptr<const jpcRes> > _tripletMm1States;
  void fillVec(boost::shared_ptr<const jpcRes>, std::vector< boost::shared_ptr<const jpcRes> >&);
  std::vector< boost::shared_ptr<const jpcRes> > extractJPCStates(std::vector< boost::shared_ptr<const jpcRes> > , std::vector< boost::shared_ptr<const jpcRes> > ) const;  
};




#endif /* _pbarpStates_H */
