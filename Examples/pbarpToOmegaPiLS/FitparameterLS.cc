#include "Examples/pbarpToOmegaPiLS/FitparameterLS.hh"

FitParameterLS::FitParameterLS()
{
}

FitParameterLS::FitParameterLS(FitParameterLS &theFitParameter)
{
  m_ParamMap = theFitParameter.m_ParamMap;
}

FitParameterLS::~FitParameterLS()
{
}

void FitParameterLS::printParamMap(std::ostream& os){

  ParameterMap::iterator it;
  
  for (it = m_ParamMap.begin(); it != m_ParamMap.end(); ++it){
    os << "\n" << it->first; 
    std::vector<double> tmpVals=it->second;
    for (unsigned int i=0; i<tmpVals.size(); ++i){
      os << "\t" << tmpVals[i];
    }
  }
  os << "\n"; 
}
