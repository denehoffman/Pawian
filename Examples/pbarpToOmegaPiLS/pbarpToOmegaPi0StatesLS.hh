#ifndef _pbarpToOmegaPi0StatesLS_H
#define _pbarpToOmegaPi0StatesLS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "PwaUtils/DataUtils.hh"


class pbarpToOmegaPi0StatesLS : public pbarpStatesLS {

public:

  /// Default Constructor  
  pbarpToOmegaPi0StatesLS(int Lmax);

  /** Destructor */
  virtual ~pbarpToOmegaPi0StatesLS();

  virtual std::vector< boost::shared_ptr<const JPCLSls> > all_JPCLSls_States() const {return _allJPCLSlsStates;};
  virtual std::vector< boost::shared_ptr<const JPCLSls> > singlet_JPCLSls_States() const {return _singlet_JPCLS_StatesDec;};
  virtual std::vector< boost::shared_ptr<const JPCLSls> > triplet0_JPCLSls_States() const {return _triplet0_JPCLS_StatesDec;};
  virtual std::vector< boost::shared_ptr<const JPCLSls> > tripletp1_JPCLSls_States() const {return _tripletp1_JPCLS_StatesDec;};
  virtual std::vector< boost::shared_ptr<const JPCLSls> > tripletm1_JPCLSls_States() const {return _tripletm1_JPCLS_StatesDec;};

  virtual void print(std::ostream& os) const;


 protected:
  virtual bool calcStates();

//   std::vector< boost::shared_ptr<const JPCLSJJls> > _allStatesDec;
  std::vector< boost::shared_ptr<const JPCLSJJls> > _lamUpUp_JPCLSJJ_StatesDec;
  std::vector< boost::shared_ptr<const JPCLSJJls> > _lamUpDown_JPCLSJJ_StatesDec;
  std::vector< boost::shared_ptr<const JPCLSJJls> > _lamDownUp_JPCLSJJ_StatesDec;
  std::vector< boost::shared_ptr<const JPCLSJJls> > _lamDownDown_JPCLSJJ_StatesDec;

  std::vector< boost::shared_ptr<const JPCLSls> > _singlet_JPCLS_StatesDec;
  std::vector< boost::shared_ptr<const JPCLSls> > _triplet0_JPCLS_StatesDec;
  std::vector< boost::shared_ptr<const JPCLSls> > _tripletp1_JPCLS_StatesDec;
  std::vector< boost::shared_ptr<const JPCLSls> > _tripletm1_JPCLS_StatesDec;

  std::vector< boost::shared_ptr<const JPCLS> > _allOmegaPiStates;
  std::vector< boost::shared_ptr<const JPCLSls> > _allJPCLSlsStates;
 private:
  jpcRes _omegaJPC;
  jpcRes _piJPC;

};




#endif /* _pbarpToOmegaPi0StatesLS_H */
