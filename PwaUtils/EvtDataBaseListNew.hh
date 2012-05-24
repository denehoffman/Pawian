#ifndef _EvtDataBaseListNew_H
#define _EvtDataBaseListNew_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <boost/unordered_set.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

typedef std::map<std::string, Vector4<double> > mapStr4Vec;
typedef std::map<std::string, map<Spin,map<Spin,map<Spin,complex<double> > > > > mapStrSpinComplex;

// typedef boost::unordered_set<std::string, Vector4<double> > mapStr4Vec;
// typedef boost::unordered_set<std::string, map<Spin,map<Spin,map<Spin,complex<double> > > > > mapStrSpinComplex;

struct EvtDataNew
  {
    mapStr4Vec FourVecsProd;
    mapStr4Vec FourVecsDec;
    mapStrSpinComplex WignerDsProd;
    mapStrSpinComplex WignerDsDec;
    double evtWeight;
  };


class EventList;

class EvtDataBaseListNew {

public:

  // create/copy/destroy:

  ///Constructor 
  EvtDataBaseListNew();



  /** Destructor */
  virtual ~EvtDataBaseListNew();

  // Getters:
  const std::vector<EvtDataNew*> getDataVecs() const {return _evtDataList;}
  const std::vector<EvtDataNew*> getMcVecs() const {return _mcDataList;}

protected:
  std::vector<EvtDataNew*> _evtDataList;
  std::vector<EvtDataNew*> _mcDataList;
  virtual void read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList)=0;  

private:

};

  
#endif /* _EvtDataBaseListNew_H */
