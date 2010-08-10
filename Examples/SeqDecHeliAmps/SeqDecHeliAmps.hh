#ifndef _SeqDecHeliAmps_H
#define _SeqDecHeliAmps_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

using namespace std;
struct jpRes
{
  Spin J;
  int P;
};
struct jpChain
{
  jpRes mother;
  jpRes daughter1;
  jpRes daughter2; 
};

struct heliAmp
{
  Spin J;
  Spin M;
  Spin lambda1;
  Spin lambda2;
  complex<double> amp;
};

class SeqDecHeliAmps {

public:

  /// Default Constructor (0x0)
  SeqDecHeliAmps(jpChain& motherChain, jpChain& daughterChain, std::vector<double>& motherMs, bool verbose);

  /** Destructor */
  virtual ~SeqDecHeliAmps();

 protected:

 private:
  void get4Vecs (std::ifstream& inStream, Vector4<double>& the4Vec);
  double EvtDecayAngle(const Vector4<double>& parent,const Vector4<double>& resonance, const Vector4<double>& daughter);
  Vector4<double> getHeliVec(const Vector4<double>& cm, const Vector4<double>& mother, const Vector4<double>& daughter);
  void calcHeliAmps(std::vector<heliAmp>& theAmps, const jpChain& theChain,  const std::vector<double>& Ms, const Vector4<double>& reference4V, const Vector4<double>& daughter1_4V, const Vector4<double>& daughter2_4V);

  std::string writeHeliAmpString(heliAmp& theHeliAmp);

  double maxAbsM(const std::vector<double>& Ms);
  void printHeliAmp(const heliAmp&);
  bool checkValidity(const jpChain& motherChain, const jpChain& daughterChain, std::vector<double>& motherMs); 
  bool _verbose;
};




#endif /* _SeqDecHeliAmps_H */
