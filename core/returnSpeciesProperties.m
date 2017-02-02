function[cp, h, s, g] = returnSpeciesProperties(T, p, Y, a, A, MW)
%	This is a function that emulates CHEMKIN, returning species properties as a function of thermodynamic state, composition, NASA polynomials (a and A) and molecular weight (MW)
	R_univ = 8.31446E7;
	T2 = power(T,2);
	T3 = power(T,3);
	T4 = power(T,4);
	Tinv = power(T,-1);
	Tlog = log(T);
	Nspecies = length(Y);
	cpvec = zeros(Nspecies, 1);
	hvec = zeros(Nspecies, 1);
	svec = zeros(Nspecies, 1);
	smix = 0;
	for i = 1:Nspecies
		if (T >= 1000) 
			cpvec(i) 	=       a(i,1) +      T.*a(i,2) +      T2.*a(i,3) +      T3.*a(i,4) +      T4.*a(i,5);
			hvec(i) 	=       a(i,1) + 0.5.*T.*a(i,2) + 1/3.*T2.*a(i,3) + 1/4.*T3.*a(i,4) + 1/5.*T4.*a(i,5) + Tinv.*a(i,6);
			svec(i) 	= Tlog.*a(i,1) +      T.*a(i,2) + 1/2.*T2.*a(i,3) + 1/3.*T3.*a(i,4) + 1/4.*T4.*a(i,5) +       a(i,7);
		else
			cpvec(i) 	=       A(i,1) +      T.*A(i,2) +      T2.*A(i,3) +      T3.*A(i,4) +      T4.*A(i,5);
			hvec(i) 	=       A(i,1) + 0.5.*T.*A(i,2) + 1/3.*T2.*A(i,3) + 1/4.*T3.*A(i,4) + 1/5.*T4.*A(i,5) + Tinv.*A(i,6);
			svec(i) 	= Tlog.*A(i,1) +      T.*A(i,2) + 1/2.*T2.*A(i,3) + 1/3.*T3.*A(i,4) + 1/4.*T4.*A(i,5) +       A(i,7);
		end
		cpvec(i) 	= cpvec(i)*R_univ/MW(i);
		hvec(i)		= hvec(i)*R_univ/MW(i)*T;
		if (Y(i) >  1E-61)
%			svec(i) 	= svec(i)*(R_univ/MW(i)) - R_univ/MW(i)*log(Y(i)*p/1E6);
			svec(i) 	= svec(i)*(R_univ/MW(i)) - R_univ/MW(i)*log(Y(i)*p/1E5);
		end
	end
	myR = R_univ./MW;
	cp 	= dot(cpvec, Y);
	h	= dot(hvec, Y);
%	      Normal Enthalpy [s0(t)] &  Mixing Rule referenced to one bar (in cgs);
	s	= dot(svec, Y)            -  smix;
	g = h - T*s;
end%returnSpeciesProperties()
