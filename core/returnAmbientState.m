function[gamma, T0, p0, Zbar] = returnAmbientState()
	global mechanism;
	if mechanism == 1
		Zbar = 0.0627964;%stoichiometric value for dodecane
		T0 = 2100.0;%
		p0 = 1E5;
		[Nspecies, species, a, A, MW, Hover] = speciesPropsC12H26();
		[gamma] = returnGamma(T0, p0, Zbar);
	elseif (mechanism == 2)
		Zbar = 0.0627964;%stoichiometric value for methane
		T0 = 2100.0;%
		p0 = 1E5;
		[Nspecies, species, a, A, MW, Hover] = speciesPropsCH4();
		[gamma] = returnGamma(T0, p0, Zbar);
	elseif (mechanism == 3)
		Zbar = 0.0627964;%stoichiometric value for hydrogen
		T0 = 2100.0;%
		p0 = 1E5;
		[Nspecies, species, a, A, MW, Hover] = speciesPropsH2();
		[gamma] = returnGamma(T0, p0, Zbar);
	elseif (mechanism == 4)
		Zbar = 0.5; %operating conditions for inerts
		T0 = 316.483;%K
		p0 = 98043.4487;%Pa
		[Nspecies, species, a, A, MW, Hover] = speciesPropsInert();
		[gamma] = returnGamma(T0, p0, Zbar);
	elseif (mechanism == 5)%Operating conditions for MW test
		Zbar = 0.0510;
		T0 = 1892;
		p0 = 2.42E5;
		[Nspecies, species, a, A, MW, Hover] = speciesPropsMW();
		[gamma] = returnGamma(T0, p0, Zbar);
	elseif (mechanism == 6)%Jeonglae algebraic psi checks
		Zbar 	= 0.5;
		T0 		= 2100;
		p0		= 1E5;
		gamma 	= 1.4;
	else
		error('mechanism not defined');
	end
end
