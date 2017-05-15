function[w, s] = charsFromInvrnts(eta, I, flaggo,SPLINES)
	[gamma] = returnAmbientState();

	if (flaggo == 0)
		M			= ppval(SPLINES(1), eta);
		Psibar		= ppval(SPLINES(4), eta);
	elseif (flaggo == 1)
		M			= MFromEtaLVG(eta, SPLINES);
		[Psibar] 	= BaseFlowFromMLVG(eta);
	else
		error('Add MfromA stuff here');	
	end

%	Pre-compute for convenience
	gm1 	= gamma - 1;
	gp1 	= gamma + 1;
	gm1o2 	= gm1/2;
	alpha 	= 1/(1+gm1o2*M*M);

%	Obtain the dimensionless perturbations at the throat and store them into param		
	alpha = 1/(1 + gm1o2*M*M);				% Common prefactor pre-computed for convenience

%	The perturbation invariants are related to the normalized perturbation quantities s = [p'/gamma *p0, u'/u, sigma'/c_p, xi] by the matrix P
	P = [  1	        1		 		 -1		  -Psibar  ;
		 gm1*alpha		alpha*gm1*M*M	alpha	alpha*Psibar;
		   0			0		  		  1			0	  ;
		   0			0		  		  0			1     ];
				
	s = P\I;

%	The characteristics w = (pi^+, pi^-, sigma, xi) are related to the normalized perturbation quantities s by the matrix R
	R = [1 	M 	0	0;
	   	 1 -M 	0 	0;
		 0	0	1	0;
		 0	0	0	1];

	w = R*s;
end
