// pbarpStatesLS class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"


class pbarpStatesLS : public AbsStates {

public:

  pbarpStatesLS();
  pbarpStatesLS(int lmax);
  virtual ~pbarpStatesLS();

  virtual std::vector< boost::shared_ptr<const jpcRes> > jpcStates() const {
    return _alljpcRes;
  }
  virtual std::vector< boost::shared_ptr<const JPCLS> > all_JPCLS_States() const {
    return _JPCLS_AllStates;
  }
  virtual std::vector< boost::shared_ptr<const JPCLS> > lamUpUp_JPCLS_States() const {
    return _lamUpUp_JPCLS_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLS> > lamUpDown_JPCLS_States() const {
    return _lamUpDown_JPCLS_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLS> > lamDownUp_JPCLS_States() const {
    return _lamDownUp_JPCLS_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLS> > lamDownDown_JPCLS_States() const {
    return _lamDownDown_JPCLS_States;
  }

  virtual std::vector< boost::shared_ptr<const JPCLSJJ> > all_JPCLSJJ_States() const {
    return _allStates;
  }
  virtual std::vector< boost::shared_ptr<const JPCLSJJ> > lamUpUp_JPCLSJJ_States() const {
    return _lamUpUp_JPCLSJJ_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLSJJ> > lamUpDown_JPCLSJJ_States() const {
    return _lamUpDown_JPCLSJJ_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLSJJ> > lamDownUp_JPCLSJJ_States() const {
    return _lamDownUp_JPCLSJJ_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLSJJ> > lamDownDown_JPCLSJJ_States() const {
    return _lamDownDown_JPCLSJJ_States;
  }

  virtual std::vector< boost::shared_ptr<const JPCLS> > singlet_JPCLS_States() const {
    return _singlet_JPCLS_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLS> > triplet0_JPCLS_States() const {
    return _triplet0_JPCLS_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLS> > tripletp1_JPCLS_States() const {
    return _tripletp1_JPCLS_States;
  }
  virtual std::vector< boost::shared_ptr<const JPCLS> > tripletm1_JPCLS_States() const {
    return _tripletm1_JPCLS_States;
  }
 
  virtual void print(std::ostream& os) const;

 protected:
  int _lmax;
  jpcRes _pbarJPC;
  jpcRes _pJPC;

  virtual bool calcStates();

  std::vector< boost::shared_ptr<const jpcRes> > _alljpcRes;
  std::vector< boost::shared_ptr<const JPCLSJJ> > _allStates;
  std::vector< boost::shared_ptr<const JPCLSJJ> > _lamUpUp_JPCLSJJ_States;
  std::vector< boost::shared_ptr<const JPCLSJJ> > _lamUpDown_JPCLSJJ_States;
  std::vector< boost::shared_ptr<const JPCLSJJ> > _lamDownUp_JPCLSJJ_States;
  std::vector< boost::shared_ptr<const JPCLSJJ> > _lamDownDown_JPCLSJJ_States;

  std::vector< boost::shared_ptr<const JPCLS> > _lamUpUp_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLS> > _lamUpDown_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLS> > _lamDownUp_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLS> > _lamDownDown_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_AllStates;
 
  std::vector< boost::shared_ptr<const JPCLS> > _singlet_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLS> > _triplet0_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLS> > _tripletp1_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLS> > _tripletm1_JPCLS_States;
};
