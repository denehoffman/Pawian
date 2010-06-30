#ifndef _OmegaPiEventList_H
#define _OmegaPiEventList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"

using OmegaPiData::OmPiEvtData;

class EventList;

class OmegaPiEventList {

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiEventList(EventList& evtListData, EventList& evtListMc);



  /** Destructor */
  virtual ~OmegaPiEventList();

  // Getters:
  const std::vector<OmPiEvtData>& getDataVecs() const {return _dataList;}
  const std::vector<OmPiEvtData>& getMcVecs() const {return _mcList;}
//   int kindOfData() const {return _kindOfData;}

protected:


private:
//   int _nOfData;
//   int _nOfMc;

  std::vector<OmPiEvtData> _dataList;
  std::vector<OmPiEvtData> _mcList;

  void read4Vecs(EventList& evtList, std::vector<OmPiEvtData>& omPiEvtList);
  
//   void get4Vecs (std::ifstream& inStream, Vector4<double>& the4Vec);

//   Tensor<complex<double> > calcSpin2Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm);
//   Tensor<complex<double> > calcSpin1Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm);

};

#endif /* _EtacToapi0EventList_H */
