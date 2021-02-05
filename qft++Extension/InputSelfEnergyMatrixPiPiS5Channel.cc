#include "qft++Extension/InputSelfEnergyMatrixPiPiS5Channel.hh"
#include "ErrLogger/ErrLogger.hh"

void InputSelfEnergyMatrixPiPiS5Channel::cacheAmps(unsigned int i, unsigned int j){
  	if(nCalled!=-99 && i<3 && j<3){
		std::string file = std::getenv("TOP_DIR");
		file = file + "/qft++Extension/SelfEnergyData/SelfEnergyMatrix_" + std::to_string(i+1) + std::to_string(j+1) + ".dat";
		inInputSelfEnergyMatrixPiPiS5Channel.open(file.c_str());

    	std::string line;
		std::string s1, s2, s3;
		while(inInputSelfEnergyMatrixPiPiS5Channel.is_open() && std::getline(inInputSelfEnergyMatrixPiPiS5Channel,line)){
		std::stringstream ss;
		ss << line;
		ss >> s1 >> s2 >> s3;
		_s.push_back(atof(s1.c_str())*atof(s1.c_str())); //table is in sqrt(s)
		_realp.push_back(atof(s2.c_str()));
		_imagp.push_back(atof(s3.c_str()));
		}
		inInputSelfEnergyMatrixPiPiS5Channel.close();
		nCalled = -99;
	}
	else if (i==3 && j==3){
		m_1=0.547862;
		m_2=0.547862;
		useAnalyticForm=true;
	}
	else if (i==4 && j==4){
		m_1=0.547862;
		m_2=0.95778;
		useAnalyticForm=true;
	}
}

complex<double> InputSelfEnergyMatrixPiPiS5Channel::interpolate(double current_s){
	complex<double>amp=0;
	if(!useAnalyticForm){

        amp = AbsComplexFunction::interpolate(current_s);

	}
	else if(useAnalyticForm){
		amp = SE4455(current_s,m_1,m_2);
	}

	return amp;
}

double InputSelfEnergyMatrixPiPiS5Channel::kaellen(double current_s, double m1, double m2){
	double out=0;
	out = (pow(m1-m2,2)-current_s)*(pow(m1+m2,2)-current_s);
	return out;
}

complex<double> InputSelfEnergyMatrixPiPiS5Channel::SE4455(double current_s, double m1, double m2){
	double par = 0;
	if(abs(m1-m2)>10e-3) par = ( (pow(m1,2)-pow(m2,2))/current_s - ( (pow(m1,2)+pow(m2,2)) / (pow(m1,2)-pow(m2,2)) ) ) * log(m1/m2) -1 ; 
	
	double nu = kaellen(current_s,m1,m2);
	double re=0;
	double im=0;
	complex<double> out=0;
	complex<double> cmplx_i(0.,1.);
	if(nu > 0){
		re = par  + sqrt(nu)/current_s * log((current_s-pow(m1,2)-pow(m2,2)+sqrt(nu))/(current_s-pow(m1,2)-pow(m2,2)-sqrt(nu)));
		im = sqrt(nu)/current_s * (-2) * PawianConstants::pi * ((current_s>pow(m1+m2,2))?1:0);
		complex<double> temp(re,im);
		out = temp;
		out*= -4*m1*m2/(16*PawianConstants::pi);
	}
	else if(nu < 0){
	  out = -4*m1*m2/(16*PawianConstants::pi) * (par + sqrt(-nu)/current_s * (atan((pow(m1,2)-pow(m2,2)+current_s)/sqrt(-nu))-atan((pow(m1,2)-pow(m2,2)-current_s)/sqrt(-nu))));
	}
	else
	{
		out = 0;
	}
	return out;
}

complex<double> InputSelfEnergyMatrixPiPiS5Channel::eval(double current_s){

    return interpolate(current_s);
}

//complex<double> InputSelfEnergyMatrixPiPiS5Channel::eval(complex<double> current_s){
//    return AbsComplexFunction::interpolate(current_s);
//    return interpolate(current_s);
//}
