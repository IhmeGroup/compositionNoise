function[cp, h, s, g, h0, gamma] = returnSpeciesProperties(T, p, Y, a, A, MW, Hover)
%	This is a function that emulates CHEMKIN, returning species properties as a function of thermodynamic state, composition, NASA polynomials (a and A) and molecular weight (MW)
	R_univ = 8.31446;
	T2 = power(T,2);
	T3 = power(T,3);
	T4 = power(T,4);
	Tinv = power(T,-1);
	Tlog = log(T);
	Nspecies = length(Y);
	cpvec = zeros(Nspecies, 1);
	hvec = zeros(Nspecies, 1);
	svec = zeros(Nspecies, 1);
	MW_bar	= 1./sum(Y./MW);
			p_bar_p_0	= p/1E5;
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
%		Following the mixing rule explained in 11.16-11.17 of Borgnakke and Sonntag
%		S = m_a*s_a + m_b*s_b
%		s_i = s_o,i + c_p*log(T/T_o) - R_i*log(p_i/p_o)
%		where the first two terms are baked into svec and the pressure adjustment is applied below
		if (Y(i) >  1E-61)%Trying to avoid the singularity associated with log(0)
			p_i_p_bar 	= Y(i)*MW_bar/MW(i);
			svec(i) 	= svec(i)*R_univ/MW(i) - R_univ/MW(i)*log(p_i_p_bar*p_bar_p_0);
		else
			svec(i)		= 0;
		end
	end
	cp 	= dot(cpvec, Y);
	h	= dot(hvec, Y);
	s	= dot(svec, Y);
	g = h - T*s;
	h0 = dot(Y,Hover./MW.*R_univ);
	gamma = cp/(cp - R_univ/MW_bar);	
end%returnSpeciesProperties()
