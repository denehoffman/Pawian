#ifndef _OmegaPiEventList_H
#define _OmegaPiEventList_H

#include <iostream>
#include <vector>

#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/text_iarchive.hpp> 

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/MATpbarpToOmegaPi/OmegaPiData.hh"

using OmegaPiData::OmPiEvtData;

class EventList;

class OmegaPiEventList {

  friend class boost::serialization::access; 

  template <typename Archive> 
  void serialize(Archive & ar, const unsigned int) {
    using boost::serialization::make_nvp;

   // unsigned _jmax;
   // unsigned _pbarmom;
   // std::vector<OmPiEvtData> _dataList;
   // std::vector<OmPiEvtData> _mcList;

    ar & BOOST_SERIALIZATION_NVP(_jmax);
    ar & BOOST_SERIALIZATION_NVP(_pbarmom);
    ar & make_nvp("dataList",_dataList);
    ar & make_nvp("mcList",_mcList);

    /* Add your own class-variables here in the following way:
      ar & BOOST_SERIALIZATION_NVP(myVar);
      or
      ar & make_nvp("myVar", myVar); // The latter form can be necessary when dealing with templates
    */
  }

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiEventList(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom);



  /** Destructor */
  virtual ~OmegaPiEventList();

  // Getters:
  const std::vector<OmPiEvtData*>& getDataVecs() const {return _dataList;}
  const std::vector<OmPiEvtData*>& getMcVecs() const {return _mcList;}
  const unsigned jMax() const {return _jmax;}
  const unsigned pbarMom() const {return _pbarmom;}

//   int kindOfData() const {return _kindOfData;}

protected:


private:

  OmegaPiEventList();

  unsigned _jmax;
  unsigned _pbarmom;
  std::vector<OmPiEvtData*> _dataList;
  std::vector<OmPiEvtData*> _mcList;

  void read4Vecs(EventList& evtList, std::vector<OmPiEvtData*>& omPiEvtList);
  void myWigner_D(const Spin &__jmax,double __alpha,double __beta,double __gamma,
	      map<serSpin,map<serSpin,map<serSpin,complex<double> > > > &__D);
  
//   void get4Vecs (std::ifstream& inStream, Vector4<double>& the4Vec);

//   Tensor<complex<double> > calcSpin2Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm);
//   Tensor<complex<double> > calcSpin1Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm);

};

#endif /* _EtacToapi0EventList_H */
