#ifndef _pbarpStates_H
#define _pbarpStates_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"


class pbarpStates : public AbsStates {

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

  std::vector< boost::shared_ptr<const jpcRes> > _singletStates;
  std::vector< boost::shared_ptr<const jpcRes> > _tripletM0States;
  std::vector< boost::shared_ptr<const jpcRes> > _tripletMp1States;
  std::vector< boost::shared_ptr<const jpcRes> > _tripletMm1States;
  void fillVec(boost::shared_ptr<const jpcRes>, std::vector< boost::shared_ptr<const jpcRes> >&);
  std::vector< boost::shared_ptr<const jpcRes> > extractJPCStates(std::vector< boost::shared_ptr<const jpcRes> > , std::vector< boost::shared_ptr<const jpcRes> > ) const;  
};




#endif /* _pbarpStates_H */
