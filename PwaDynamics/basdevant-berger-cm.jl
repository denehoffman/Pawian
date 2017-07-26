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
			tmp = (((1.0/pi)*fR2^2*imag(Sigma(a2,x[2], m12, m22)))/abs2(ds(x[2],m12, m22, fR2, mR2, a2)))*cm(sreal, simag, sqrt(sprime(a2, x[2])), mR1)*(1.0/((1.0-x[2])^2))
		elseif fR2 == 0.0
			tmp = (((1.0/pi)*fR1^2*imag(Sigma(a1,x[1], m11, m21)))/abs2(ds(x[1],m11, m21, fR1, mR1, a1)))*cm(sreal, simag, sqrt(sprime(a1, x[1])), mR2)*(1.0/((1.0-x[1])^2))
		else
			tmp = ((1.0/(pi^2))*fR1^2*imag(Sigma(a1,x[1], m11, m21)))/abs2(ds(x[1],m11, m21, fR1, mR1, a1))*cm(sreal, simag, sqrt(sprime(a1, x[1])), sqrt(sprime(a2, x[2])))*(1.0/((1.0-x[1])^2))*(1.0/((1.0-x[2])^2))*fR2^2*imag(Sigma(a2,x[2], m12, m22))/abs2(ds(x[2],m12, m22, fR2, mR2, a2))
		end
		f[1], f[2] = reim(tmp)
	end
	if fR1 == 0.0 && fR2!=0.0
		result = cuhre(integrand,1,2)
	elseif fR2 == 0.0 && fR1!=0.0
		result = cuhre(integrand,1,2)
	elseif fR2 == 0.0 && fR1==0.0
		result = cm(sreal, simag, mR1, mR2)
	else
		result = cuhre(integrand,2,2)
	end
	
	resultCmplx = Float64[ real(complex(result[1]...)), imag(complex(result[1]...)) ]
	
	return resultCmplx
end
end

using PhaseSpace
for i = 0:400
	for j = 0:0
		iscale = i/100.0
		jscale = j/100.0
		#result = bbcm(iscale,jscale,0.14,0.14,0.14,0.14,0.758287544,0.442718872,0.14,0.0,0.0)
		result = bbcm(iscale,jscale,0.14,0.14,0.14,0.14,0.77,0.0,0.14,0.0,0.0)
		betrag=abs(complex(result[1],result[2]))
		println(sqrt(iscale), " ",sqrt(jscale), " ",  betrag)
	end
end


