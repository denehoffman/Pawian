#ifndef MIGENERIC_HH_
#define MIGENERIC_HH_

#include <vector>

class MIGeneric
{
public:
  MIGeneric()
	  {
	  }

  virtual ~MIGeneric()
	{ /* nothing */	}

  // Evaluate a given parameter set
  virtual double evaluator() = 0;
  // Retrieve a clone of a GenericInterface object-derivative
  // Needed so that we can clone objects through a base-pointer
  // without information on whether we are dealing with PAWIAN or GPGPU.
  virtual boost::shared_ptr<MIGeneric> clone() = 0;

  // some access code for adding, modifying and retrieving parameters
  virtual void setStartPar(int num, double *par, double* min, double* max, double* err) = 0;
  virtual const double getPar(const unsigned int num) const {
    if(num>=parValue_.size()) return 0;
    return parValue_[num];
  };
  virtual const double getMin(const unsigned int num) const {
    if(num>=minValue_.size()) return 0;
    return minValue_[num];
  };
  virtual const double getMax(const unsigned int num) const {
    if(num>=maxValue_.size()) return 0;
    return maxValue_[num];
  };
  virtual const double getErr(const unsigned int num) const {
    if(num>=parError_.size()) return 0;
    return parError_[num];
  };
  virtual const unsigned int getNumPars() const {
    return parValue_.size();
  };

protected:
  std::vector<double> parValue_;
  std::vector<double> minValue_;
  std::vector<double> maxValue_;
  std::vector<double> parError_;
};

#endif
