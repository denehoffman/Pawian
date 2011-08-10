#ifndef MIBASE_HH_
#define MIBASE_HH_

class MIBase
{

public:

  MIBase()
	  {
	  }

  virtual ~MIBase()
	{ /* nothing */	}

  virtual const double exec(unsigned int Npar, double* par,  double* min, double* max, double* err) =0;
 
};

#endif
