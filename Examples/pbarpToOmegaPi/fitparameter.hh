//
// C++ Interface: fitparameter
//
// Description: 
//
//
// Author: Denis Pavlina <denis@pc15>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef FITPARAMETER_H
#define FITPARAMETER_H

#include <vector>
#include <map>
#include <string>

using namespace std;

//Typedef for user supplied minuit parameter map
typedef map<string, vector<double> > ParameterMap;

/**
General named fit parameter class.

	@author Denis Pavlina <denis@pc15>
*/
class FitParameter
{
  public:
    FitParameter();
    
    FitParameter(FitParameter &theFitParameter);

    ~FitParameter();
    
    inline ParameterMap& getParamMap()
    {
      return m_ParamMap;
    }

  protected:
    ParameterMap m_ParamMap;
};

#endif
