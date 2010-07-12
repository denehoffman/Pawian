#ifndef EVENT_HH
#define EVENT_HH

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
  ~Event();

  void addParticle(double e, double px, double py, double pz);
  
  Vector4<float>* p4(int i);
  float* pid(int i);

  int size();
  
private:
  std::vector<EvtPartData*> particles;
  double weight;
};

#endif
