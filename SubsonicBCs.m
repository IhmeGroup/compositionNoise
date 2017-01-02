function[Res] = SubsonicBCs(I_a, I_b)
%	Boundary conditions for subsonic flows
%
%	Unpack the parameter vector
	global param;
	M_a 	= param(1);
	M_b 	= param(2);
	M_c 	= param(3);
	gamma 	= param(4);
%	Omega	= param(5);
	Tbar	= param(6);
	pbar	= param(7);
	Zbar	= param(8);
	w_p_a	= param(9);
%	w_m_a	= param(10);
	w_s_a 	= param(11);
	w_z_a	= param(12);
%	w_p_b	= param(13);
	w_m_b	= param(14);
%	w_s_b 	= param(15);
%	w_z_b	= param(16);

%	Pre-compute for convenience
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	alpha_a = 1/(1+gm1o2*M_a*M_a);
	alpha_b = 1/(1+gm1o2*M_b*M_b);

	T_a = (1 + gm1o2*M_a*M_a)^(-1)*Tbar;
	T_b = (1 + gm1o2*M_b*M_b)^(-1)*Tbar;
	p_a = (1 + gm1o2*M_a*M_a)^(-gamma/gm1)*pbar;
	p_b = (1 + gm1o2*M_b*M_b)^(-gamma/gm1)*pbar;
	Psi_a = returnPsi(T_a, p_a, Zbar);
	Psi_b = returnPsi(T_b, p_b, Zbar);


	P_a = [ 	1 					   1 			  -1			-Psi_a   ;
			gm1*alpha_a 	gm1*M_a*M_a*alpha_a 	alpha_a		alpha_a*Psi_a;
				0 					   0 			   1			   0     ;
				0					   0			   0			   1     ];

	P_b = [ 	1 					   1 			  -1			-Psi_b   ;
			gm1*alpha_b 	gm1*M_b*M_b*alpha_b 	alpha_b		alpha_b*Psi_b;
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
