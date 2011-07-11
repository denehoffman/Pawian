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
#ifndef FITPARAMETERLS_H
#define FITPARAMETERLS_H

#include <iostream>
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
class FitParameterLS
{
  public:
    FitParameterLS();
    
    FitParameterLS(FitParameterLS &theFitParameter);

    ~FitParameterLS();
    
    inline ParameterMap& getParamMap()
    {
      return m_ParamMap;
    }

   void printParamMap(std::ostream& os);

  protected:
    ParameterMap m_ParamMap;
};

#endif
