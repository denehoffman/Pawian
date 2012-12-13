#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class IsobarDecay;

class pbarpReaction {

public:
  pbarpReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, int lmax=3);

  virtual ~pbarpReaction();

  virtual void print(std::ostream& os) const;
  

protected:

private:
  int _lmax;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpJPCs;
  std::vector< boost::shared_ptr<IsobarDecay> > _prodDecs;

  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarDecay> >, pawian::Collection::SharedPtrLess > _pbarpSingletDecMap;
  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarDecay> >, pawian::Collection::SharedPtrLess > _pbarpTriplet0DecMap;
  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarDecay> >, pawian::Collection::SharedPtrLess > _pbarpTripletp1DecMap;
  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarDecay> >, pawian::Collection::SharedPtrLess > _pbarpTripletm1DecMap;

  void  fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<IsobarDecay> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarDecay> >, pawian::Collection::SharedPtrLess > toFill);
};
