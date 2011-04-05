//
// C++ Interface: spindensity
//
// Description: 
//
//
// Author: Denis Pavlina <denis@pc15>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SPINDENSITY_H
#define SPINDENSITY_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"
#include "PwaUtils/DataUtils.hh"

// using namespace std;


class OmegaPiEventList;
class pbarpToOmegaPi0States;

/**
Class for calculation of spin density matrix

	@author Denis Pavlina <denis@pc15>
*/
class spinDensity
{
public:
    static complex<double> calcSpinDensityIncoherent(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal);
    static complex<double> calcSpinDensityOmegaFrame(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal);
    static complex<double> calcIncoherentAmp(Spin Minit, Spin lamomega, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData);


    spinDensity();

    ~spinDensity();

};

#endif
