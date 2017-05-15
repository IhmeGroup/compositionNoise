function[d_I_d_eta] = DuranMoreauODE(eta, I, Omega, flaggo, SPLINES)	
%	This is a function that returns the rhs of the spatially evolving ODE as a function of axial coordinate for nozzles with a linear velocity profile

	[gamma, T0, p0, Zbar] = returnAmbientState();
%	Precompute for convenience
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;

%	Unpack the spline vector
	if (flaggo == 0)
		M 		= ppval(SPLINES(1), eta);
		M2 		= M*M;
		Psibar 	= ppval(SPLINES(4), eta);
		ubar 	= ppval(SPLINES(5), eta);
	elseif(flaggo == 1)
		M				= MFromEtaLVG(eta, SPLINES);
		M2				= M*M;
		[Psibar, ubar] 	= BaseFlowFromMLVG(M);
	else
		error('Add MfromA stuff here');	
	end

	A = -2*pi*sqrt(-1)*Omega/(ubar*(M*M-1))*[ 			M2			-(1+gm1o2*M2)/gm1				gamma/gm1					Psibar*gamma/gm1;
											 -gm1*M2/(1+gm1o2*M2)			M2				-(gm1*M2+1)/(1+gm1o2*M2)	-(1+gm1*M2)/(1+gm1o2*M2)*Psibar;
														0					0							M2-1							0;
														0					0							0								(M2-1)];

	d_I_d_eta = A*I;


end
