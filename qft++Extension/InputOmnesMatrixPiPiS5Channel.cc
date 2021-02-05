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
    _s.push_back(atof(s1.c_str())*atof(s1.c_str())); //table is in sqrt(s)
	_realp.push_back(atof(s2.c_str()));
	_imagp.push_back(atof(s3.c_str()));
	}
	inInputOmnesMatrixPiPiS5Channel.close();
	nCalled = -99;
  }
}
