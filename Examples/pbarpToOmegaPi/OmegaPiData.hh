#ifndef _OmegaPiData_H
#define _OmegaPiData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

// struct fitParamVal
// {
//   double aMass;
//   double aWidth;
//   double cont0spin;
//   double cont1spin;
//   double cont2spin;
// };

struct OmPiEvtData
{
  Vector4<float> cm_4Vec;
  Vector4<float> omegaHeliCm4Vec;
  Vector4<float> pi0RecHeliCm4Vec;
  Vector4<float> pi0HeliOmega4Vec;
//   Tensor<complex<double> > spin2_0_amp;
//   Tensor<complex<double> > spin2_1_amp;
//   Tensor<complex<double> > spin1_0_amp;
//   Tensor<complex<double> > spin1_1_amp;
  
};

#endif /* _etacToapi0Data_H */
