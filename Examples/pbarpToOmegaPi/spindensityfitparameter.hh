//
// C++ Interface: spindensityfitparameter
//
// Description: 
//
//
// Author: Denis Pavlina <denis@pc15>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SPINDENSITYFITPARAMETER_H
#define SPINDENSITYFITPARAMETER_H

#include "Examples/pbarpToOmegaPi/minuitstartparam.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"
//#include "PwaUtils/DataUtils.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "ErrLogger/ErrLogger.hh"
/**
Class for the fit parameters used to calculate spin density matrix

	@author Denis Pavlina <denis@pc15>
*/
class SpinDensityFitParameter : public minuitStartParam
{
public:
    SpinDensityFitParameter();

    ~SpinDensityFitParameter();
    
    inline void FillFitParamVal(std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > &fitParmS,
                                const std::vector< boost::shared_ptr<const JPCLS> > &theJPCLS,
                                const string &theSuffix
                               )
    { 
      std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
      for ( itJPCLS=theJPCLS.begin(); itJPCLS!=theJPCLS.end(); ++itJPCLS)
      {
        string strName = (*itJPCLS)->name()+theSuffix;
        ParameterMap::iterator it;

        it = m_ParamMap.find(strName);
        if (it != m_ParamMap.end()) 
        {
          double dMag = it->second[0];
          double dPhi = it->second[1];
          std::pair <double,double> tmpParameter=make_pair(dMag,dPhi);
          fitParmS[(*itJPCLS)]=tmpParameter; 
        }
        else
        {
          Warning <<"parameter " << strName <<" not available in fitVal input file!!!!"  << endmsg;
        }
      }
    }
    
    void getFitParamVal(/*OUT*/ OmegaPiData::fitParamVal& theFitParamVal, 
                        /*IN*/ boost::shared_ptr<const pbarpToOmegaPi0States> thePi0States
                       );
};

#endif
