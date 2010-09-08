#ifndef _OmegaPiEventList_H
#define _OmegaPiEventList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/MATpbarpToOmegaPi/OmegaPiData.hh"

using OmegaPiData::OmPiEvtData;

class EventList;

class OmegaPiEventList {

public:

  static boost::shared_ptr<OmegaPiEventList> getList();
  static boost::shared_ptr<OmegaPiEventList> getList(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom);
  // create/copy/destroy:
 
  //bool initList(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom);

  /** Destructor */
  virtual ~OmegaPiEventList();

  // Getters:
  const std::vector<OmPiEvtData*>& getDataVecs() const {return _dataList;}
  const std::vector<OmPiEvtData*>& getMcVecs() const {return _mcList;}
  const unsigned jMax() const {return _jmax;}
  const unsigned pbarMom() const {return _pbarmom;}

//   int kindOfData() const {return _kindOfData;}

protected:

  OmegaPiEventList(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom);

private:

  static boost::shared_ptr<OmegaPiEventList> theList;

  //bool initialized;
  unsigned _jmax;
  unsigned _pbarmom;
  std::vector<OmPiEvtData*> _dataList;
  std::vector<OmPiEvtData*> _mcList;

  void read4Vecs(EventList& evtList, std::vector<OmPiEvtData*>& omPiEvtList);
  
//   void get4Vecs (std::ifstream& inStream, Vector4<double>& the4Vec);

//   Tensor<complex<double> > calcSpin2Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm);
//   Tensor<complex<double> > calcSpin1Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm);

};

#endif /* _EtacToapi0EventList_H */
