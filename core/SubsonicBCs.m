function[Res] = SubsonicBCs(I_a, I_b, w_l, w_r, eta_l, eta_r, flaggo, SPLINES)
%	Boundary conditions for subsonic flows
%
%	Unpack the parameter vector
	[gamma] = returnAmbientState();
	w_p_a	= w_l(1);
	w_s_a 	= w_l(3);
	w_z_a	= w_l(4);
	w_m_b	= w_r(2);

	if (flaggo == 0)
		M_a			= ppval(SPLINES(1), eta_l);
		Psibar_a	= ppval(SPLINES(4), eta_l);
		M_b			= ppval(SPLINES(1), eta_r);
		Psibar_b	= ppval(SPLINES(4), eta_r);
	elseif (flaggo == 1)
		[M_a] = MFromEtaLVG(eta_l, SPLINES);
		[Psibar_a] = BaseFlowFromMLVG(M_a);
		[M_b] = MFromEtaLVG(eta_r, SPLINES);
		[Psibar_b] = BaseFlowFromMLVG(M_b);
	else
		error('Add MfromA stuff here');	
	end

%	Pre-compute for convenience
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	alpha_a = 1/(1+gm1o2*M_a*M_a);
	alpha_b = 1/(1+gm1o2*M_b*M_b);

	P_a = [ 	1 					   1 			  -1			-Psibar_a   ;
			gm1*alpha_a 	gm1*M_a*M_a*alpha_a 	alpha_a		alpha_a*Psibar_a;
				0 					   0 			   1			   0     ;
				0					   0			   0			   1     ];

	P_b = [ 	1 					   1 			  -1			-Psibar_b   ;
			gm1*alpha_b 	gm1*M_b*M_b*alpha_b 	alpha_b		alpha_b*Psibar_b;
				0 					   0 			   1			   0     ;
				0					   0			   0			   1     ];

	R_a = [	1 	M_a 	0	0;
		   	1  -M_a 	0	0;
			0	0		1	0;
			0	0		0	1];

	R_b = [	1 	M_b 	0	0;
		   	1  -M_b 	0	0;
			0	0		1	0;
			0	0		0	1];

	s_a = P_a\I_a;

	s_b = P_b\I_b;

	w_a = R_a*s_a;

	w_b = R_b*s_b;

	w_m_b_num = w_b(2);
	w_p_a_num = w_a(1);
	w_s_a_num = w_a(3);
	w_z_a_num = w_a(4);

	Res(1) = (w_m_b_num - w_m_b);%	(u-c) = prescribed at the outlet
	Res(2) = (w_p_a_num - w_p_a);%	(u+c) = prescribed at the inlet
	Res(3) = (w_s_a_num - w_s_a);%	(s)	  = prescirbed at the inlet
	Res(4) = (w_z_a_num - w_z_a);%	(z)	  = prescribed at the inlet
end
