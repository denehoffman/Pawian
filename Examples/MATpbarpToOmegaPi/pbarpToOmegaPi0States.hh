#ifndef _pbarpToOmegaPi0States_H
#define _pbarpToOmegaPi0States_H

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

class pbarpStates;

class pbarpToOmegaPi0States : public AbsStates {

  friend class boost::serialization::access; 

  template <typename Archive> 
  void serialize(Archive & ar, const unsigned int) {
    using boost::serialization::make_nvp;

    //ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GParameterSet);

    /*jpcRes _omegaJPC;
    jpcRes _piJPC;

    boost::shared_ptr<pbarpStates> _pbarpStatesAll;
    std::vector< boost::shared_ptr<const jpcRes> > _pbarpSingletToOmegaPi;
    std::vector< boost::shared_ptr<const jpcRes> > _pbarpTripletM0ToOmegaPi;
    std::vector< boost::shared_ptr<const jpcRes> > _pbarpTripletM1ToOmegaPi;
    std::vector< boost::shared_ptr<const JPCLS> >  _JPCLSomegaProdSinglet;
    std::vector< boost::shared_ptr<const JPCLS> > _JPCLSomegaProdTripletM0;
    std::vector< boost::shared_ptr<const JPCLS> > _JPCLSomegaProdTripletM1;*/

    ar & BOOST_SERIALIZATION_NVP(_omegaJPC);
    ar & BOOST_SERIALIZATION_NVP(_piJPC);
    ar & BOOST_SERIALIZATION_NVP(_pbarpStatesAll);
    ar & BOOST_SERIALIZATION_NVP(_pbarpSingletToOmegaPi);
    ar & BOOST_SERIALIZATION_NVP(_pbarpTripletM0ToOmegaPi);
    ar & BOOST_SERIALIZATION_NVP(_pbarpTripletM1ToOmegaPi);
    ar & BOOST_SERIALIZATION_NVP(_JPCLSomegaProdSinglet);
    ar & BOOST_SERIALIZATION_NVP(_JPCLSomegaProdTripletM0);
    ar & BOOST_SERIALIZATION_NVP(_JPCLSomegaProdTripletM1);

    /* Add your own class-variables here in the following way:
      ar & BOOST_SERIALIZATION_NVP(myVar);
      or
      ar & make_nvp("myVar", myVar); // The latter form can be necessary when dealing with templates
    */
  }

public:

  /// Default Constructor
  pbarpToOmegaPi0States(boost::shared_ptr<pbarpStates> pbarpStates);

  /** Destructor */
  virtual ~pbarpToOmegaPi0States();

  virtual std::vector< boost::shared_ptr<const jpcRes> > pbarpSingletStates() const {return _pbarpSingletToOmegaPi;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > pbarpTripletM0States() const {return _pbarpTripletM0ToOmegaPi;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > pbarpTripletM1States() const {return _pbarpTripletM1ToOmegaPi;};
  virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsSinglet() const {return _JPCLSomegaProdSinglet;};
  virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsTriplet0() const {return _JPCLSomegaProdTripletM0;};
  virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsTriplet1() const {return _JPCLSomegaProdTripletM1;};

  virtual void print(std::ostream& os) const;


 protected:
  virtual bool calcJPCs();
  vector<LS> myGetValidLS(const serSpin &__j,int __parity,const serSpin &__s1,int __p1,
		      const serSpin &__s2,int __p2);

 private:

  pbarpToOmegaPi0States();

  jpcRes _omegaJPC;
  jpcRes _piJPC;

  boost::shared_ptr<pbarpStates> _pbarpStatesAll;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpSingletToOmegaPi;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpTripletM0ToOmegaPi;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpTripletM1ToOmegaPi;
  std::vector< boost::shared_ptr<const JPCLS> >  _JPCLSomegaProdSinglet;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSomegaProdTripletM0;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSomegaProdTripletM1;
  std::vector< boost::shared_ptr<const JPCLS> > extractJPCLSStates(std::vector< boost::shared_ptr<const jpcRes> > , std::vector< boost::shared_ptr<const JPCLS> > ) const;
};




#endif /* _pbarpToOmegaPi0States_H */
