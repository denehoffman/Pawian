// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>

// Boost header files go here
#include <boost/lexical_cast.hpp>

//#include "ErrLogger/ErrLogger.hh"

// Minimizer Interface header files go here
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIMinuit.hh"
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIGeneva.hh"

// The toy-data to fit to
#include "Examples/Tutorial/FitIF/PolyFit.hh"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){

  double p0, p1, p2, p3, sigma_smear;
  p0 = -10. ; p1 = 10.; p2 = 1. ; p3 = -0.01; sigma_smear = 3;

  // Generate data distribution
  boost::shared_ptr<PolyFit> myFit(new PolyFit(p0, p1, p2, p3, sigma_smear));

  //--------------------------Minimizer IF --------------------------------------------------------

  double val[4], min[4], max[4], err[4];
  val[0] = -11; max[0] = 0; min[0] = -20; err[0] = 3;
  val[1] = 9.8; max[1] = 15; min[1] = 5; err[1] = 2;
  val[2] = 1.1; max[2] = 1.5; min[2] = 0.5; err[2] = 0.3;
  val[3] = -0.008; max[3] = 0.; min[3] = -0.02; err[3] = 0.005; 

  MIGeneva genevaFit(myFit);  
  double genResult = genevaFit.exec(4, val,  min, max, err); 

  std::cout << "final Geneva par :\t" << genResult << std::endl;

  std::cout << "final a:\t" << val[0] << " +- " << err[0] << std::endl;
  std::cout << "final b:\t" << val[1] << " +- " << err[1] << std::endl; 
  std::cout << "final c:\t" << val[2] << " +- " << err[2] << std::endl;
  std::cout << "final d:\t" << val[3] << " +- " << err[3] << std::endl; 

  //myFit->drawGraph(resultPar.at(0),resultPar.at(1),resultPar.at(2),resultPar.at(3));

  std::cout << "Done Geneva ..." << std::endl << std::endl;

  //val[0] = -11; max[0] = 0; min[0] = -20; err[0] = 3;
  //val[1] = 9.8; max[1] = 15; min[1] = 5; err[1] = 2;
  //val[2] = 1.1; max[2] = 1.5; min[2] = 0.5; err[2] = 0.3;
  //val[3] = -0.008; max[3] = 0.; min[3] = -0.02; err[3] = 0.005; 

  MIMinuit minuitFit(myFit);  
  double minResult = minuitFit.exec(4, val,  min, max, err); 

  std::cout << "final Minuit par :\t" << minResult << std::endl;

  std::cout << "final a:\t" << val[0] << " +- " << err[0] << std::endl;
  std::cout << "final b:\t" << val[1] << " +- " << err[1] << std::endl; 
  std::cout << "final c:\t" << val[2] << " +- " << err[2] << std::endl;
  std::cout << "final d:\t" << val[3] << " +- " << err[3] << std::endl; 

  std::cout<< "Done Minuit ..." << std::endl << std::endl;

   myFit->drawGraph(val[0],val[1],val[2],val[3]);

  return 0;
}
