#include "qft++Extension/InputOmnesMatrixPiPiS5Channel.hh"
#include "ErrLogger/ErrLogger.hh"

void InputOmnesMatrixPiPiS5Channel::cacheAmps(unsigned int i, unsigned int j){
  if(nCalled!=-99){
	std::string file = std::getenv("TOP_DIR");
	file = file + "/qft++Extension/OmnesData/OmnesMatrix_" + std::to_string(i+1) + std::to_string(j+1) + ".dat";
	inInputOmnesMatrixPiPiS5Channel.open(file.c_str());

  std::string line;
	std::string s1, s2, s3;
	while(inInputOmnesMatrixPiPiS5Channel.is_open() && std::getline(inInputOmnesMatrixPiPiS5Channel,line)){
	std::stringstream ss;
	ss << line;
	ss >> s1 >> s2 >> s3;
	s.push_back(atof(s1.c_str())*atof(s1.c_str())); //table is in sqrt(s)
	realp.push_back(atof(s2.c_str()));
	imagp.push_back(atof(s3.c_str()));
	}
	inInputOmnesMatrixPiPiS5Channel.close();
	nCalled = -99;
  }
}

complex<double> InputOmnesMatrixPiPiS5Channel::interpolate(double current_s){
	std::vector<double>::iterator it;
	it = std::upper_bound(s.begin(), s.end(), current_s);
	int low = (it-s.begin())-1;
	int up = (it-s.begin());	
	if(up == (int)s.size()){
		low -= 1;
		up -= 1;
	}	

	// get s, real and imag for lower and upper tabulated values
	double slow = s.at(low); double sup = s.at(up);
	double reallow = realp.at(low); double realup = realp.at(up);
	double imaglow = imagp.at(low); double imagup = imagp.at(up);
	
	// linear interpolation to actual s
	double m_real = (realup - reallow)/(sup - slow);
	double b_real = (realup + reallow-m_real*(sup + slow))/2.;
	double actual_real = m_real * current_s + b_real;

	double m_imag = (imagup - imaglow)/(sup - slow);
	double b_imag = (imagup + imaglow-m_imag*(sup + slow))/2.;
	double actual_imag = m_imag * current_s + b_imag;

	complex<double>amp = complex<double>(actual_real,actual_imag);
	return amp;
}


complex<double> InputOmnesMatrixPiPiS5Channel::interpolate(complex<double>current_s){
  return complex<double>(0,0);
}
