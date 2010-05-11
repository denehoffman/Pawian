#ifndef _etacToa1320pi0Data_H
#define _etacToa1320pi0Data_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

struct fitParamVal
{
  double a1320Mass;
  double a1320Width;
  double cont0spin;
  double cont1spin;
  double cont2spin;
};

struct evt4Vec
{
  Vector4<double> pi0_0_4Vec;
  Vector4<double> pi0_1_4Vec;
  Vector4<double> eta_4Vec;
  Vector4<double> a2_0_4Vec;
  Vector4<double> a2_1_4Vec;
  Vector4<double> cm_4Vec;
  Tensor<complex<double> > spin2_0_amp;
  Tensor<complex<double> > spin2_1_amp;
  Tensor<complex<double> > spin1_0_amp;
  Tensor<complex<double> > spin1_1_amp;  
};

#endif /* _etacToa1320pi0Data_H */
