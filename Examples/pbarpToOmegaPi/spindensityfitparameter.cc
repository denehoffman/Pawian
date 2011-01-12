#include "Examples/pbarpToOmegaPi/spindensityfitparameter.hh"

SpinDensityFitParameter::SpinDensityFitParameter()
 : minuitStartParam()
{
}


SpinDensityFitParameter::~SpinDensityFitParameter()
{
}





/*!
    \fn SpinDensityFitParameter::getFitParamVal()
 */
void SpinDensityFitParameter::getFitParamVal(/*OUT*/ OmegaPiData::fitParamVal& theFitParamVal, 
                                             /*IN*/ boost::shared_ptr<const pbarpToOmegaPi0States> thePi0States
                                            )
{
  FillFitParamVal(theFitParamVal.omegaProdSinglet,thePi0States->jpclsSinglet(),"S");
  FillFitParamVal(theFitParamVal.omegaProdTriplet0,thePi0States->jpclsTriplet0(),"T0");
  FillFitParamVal(theFitParamVal.omegaProdTriplet1,thePi0States->jpclsTriplet1(),"T1");
}
