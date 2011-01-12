#include "Examples/pbarpToOmegaPi/fitparameter.hh"

FitParameter::FitParameter()
{
}

FitParameter::FitParameter(FitParameter &theFitParameter)
{
  m_ParamMap = theFitParameter.m_ParamMap;
}

FitParameter::~FitParameter()
{
}


