#pragma once

#include "qft++/topincludes/tensor.hh"
#include "qft++/tensor/Vector4.hh"
#include "ErrLogger/ErrLogger.hh"
#include <vector>

enum pidType {proton, pion, kaon, muon, electron, photon, resonance, maxPid};

struct EvtPartData {
  Vector4<float>* vector4;
  float pidVector[maxPid];
  EvtPartData(double e, double px, double py, double pz) {
    vector4 = new Vector4<float>(e,px,py,pz);
  }
  ~EvtPartData() {
    delete vector4;
  }
};

class Event
{
public:
  Event();
  Event(int evtNumber);
  ~Event();

  void addParticle(double e, double px, double py, double pz);
  void addWeight(double);

  Vector4<float>* p4(unsigned int i);
  float* pid(unsigned int i);

  int size();
  double Weight(){ return weight; }
  int eventNo() const {return evtNo; }
  
  bool operator<(const Event& compare) const;

private:
  std::vector<EvtPartData*> particles;
  double weight;
  int evtNo;
};
