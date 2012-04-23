#include <utility>


struct dynamicModelParams{
  enum enumDynamicModel {BreitWigner, BreitWignerBlattW, Flatte, MassIndependent};
  enumDynamicModel dynamicModel;
  double mass;
  double width;
  double gFactor1;
  double gFactor2;
  pair<double, double> massPair1;
  pair<double, double> massPair2;
  complex<double> value;
  
};
