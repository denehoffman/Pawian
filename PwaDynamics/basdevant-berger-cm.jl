#************************************************************************#
#									 #
#  Copyright 2017 Markus Kuhlmann (mkuhlmann@ep1.rub.de)		 #
#  	      	   - Ruhr-Universität Bochum 				 #
#									 #
#  This file is part of Pawian.					         #
#									 #
#  Pawian is free software: you can redistribute it and/or modify	 #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or 	 #
#  (at your option) any later version.	 	      	  	   	 #
#									 #
#  Pawian is distributed in the hope that it will be useful,		 #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of	 #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	 #
#  GNU General Public License for more details.	      		         #
#									 #
#  You should have received a copy of the GNU General Public License     #
#  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	 #
#									 #
#************************************************************************#
module PhaseSpace
using Cuba
export cm, bbcm
function cm(real, imag, m, µ)
		if real == 0.0 && imag == 0.0
    			return 0.0
		end
		s=complex(real,imag)
		function A(s,m,µ)
			return (-1.0/s)*sqrt((m+µ)^2-s)*sqrt((m-µ)^2-s)
		end
		function B(s,m,µ)
			return log((sqrt((m+µ)^2-s)+sqrt((m-µ)^2-s))/(2.0*sqrt(m*µ)))
		end
		function C(s,m,µ)
			if m==µ
				return -1.0
			else
				return ((m^2-µ^2)/(2.0*s))*log(m/µ)-((m^2+µ^2)/(2.0*(m^2-µ^2)))*log(m/µ)-0.5
			end
		end
		
		return (-2.0/pi)*(A(s,m,µ)*B(s,m,µ)+C(s,m,µ))
	end
function bbcm(sreal,simag,m11,m21,m12,m22,mR1,fR1,mR2,fR2,epsilon)
	a1 = (m11+m21)^2
	a2 = (m12+m22)^2
	function scaler(u)
		return (1.0/((1.0-u)^2))
	end
	function sprime(a, t)
		return (a+(t/(1.0-t)))
	end
	function ds(t, m1, m2, fR, mR, a)
		return (sprime(a,t)-mR^2+(fR^2*(sprime(a,t)-a)*cm(sprime(a,t), epsilon ,m1,m2)))
	end
	function Sigma(a, t, m1, m2)
		return (sprime(a, t)-a)*cm(sprime(a, t), epsilon,m1,m2)
	end
	function integrand(x,f)
		if fR1 == 0.0
			if abs(x[0]-1.0) < eps(0.0)
				return 0.0
			end
			tmp = (((1.0/pi)*fR2^2*imag(Sigma(a2,x[2], m12, m22)))/abs2(ds(x[2],m12, m22, fR2, mR2, a2)))*cm(sreal, simag, sqrt(sprime(a2, x[2])), mR1)*scaler(x[2])
		elseif fR2 == 0.0
			if abs(x[1]-1.0) < eps(0.0)
				return 0.0
			end
			tmp = (((1.0/pi)*fR1^2*imag(Sigma(a1,x[1], m11, m21)))/abs2(ds(x[1],m11, m21, fR1, mR1, a1)))*cm(sreal, simag, sqrt(sprime(a1, x[1])), mR2)*scaler(x[1])
			
		else
			if abs(x[1]-1.0) < eps(0.0) || abs(x[2]-1.0) < eps(0.0)
				return 0.0
			end
			tmp = ((1.0/(pi^2))*fR1^2*imag(Sigma(a1,x[1], m11, m21)))/abs2(ds(x[1],m11, m21, fR1, mR1, a1))*cm(sreal, simag, sqrt(sprime(a1, x[1])), sqrt(sprime(a2, x[2])))*scaler(x[1])*scaler(x[2])*fR2^2*imag(Sigma(a2,x[2], m12, m22))/abs2(ds(x[2],m12, m22, fR2, mR2, a2))
		end
		f[1], f[2] = reim(tmp)
	end
	if fR1 == 0.0 && fR2 != 0.0
		result = cuhre(integrand,1,2)
	elseif fR2 == 0.0 && fR1!= 0.0
		result = cuhre(integrand,1,2)
	elseif fR2 == 0.0 && fR1 == 0.0
		result = cm(sreal, simag, mR1, mR2)
	else
		result = cuhre(integrand,2,2)
	end
	resultCmplx = Float64[ real(complex(result[1]...)), imag(complex(result[1]...)) ]
	return resultCmplx
end
end
