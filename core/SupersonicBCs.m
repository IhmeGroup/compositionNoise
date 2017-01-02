function[Res] = SupersonicBCs(I_a, I_b)
%	Compute the boundary condition residual for the supersonic side of a choked flow
%
%	Unpack the param vector
	global param;
	M_a 	= param(1);%1
	M_b 	= param(2);%2
	M_c 	= param(3);%3
	gamma 	= param(4);%4
%	Omega	= param(5);%5
	Tbar 	= param(6);%6
	pbar	= param(7);
	Zbar	= param(8);
	w_p_a	= param(9);
	w_m_a	= param(10);
	w_s_a 	= param(11);
	w_z_a 	= param(12);
%	w_p_b	= param(13);
%	w_m_b	= param(14);
%	w_s_b 	= param(15);
%	w_z_b	= param(16);

	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	alpha_a = 1/(1+gm1o2*M_a*M_a);
	T_a = (1 + gm1o2*M_a*M_a).^(-1)*Tbar;
	p_a = (1 + gm1o2*M_a*M_a).^(-gamma/gm1)*pbar;
	Psi_a = returnPsi(T_a, p_a, Zbar);

	P_a = [		1			1				  -1	  -Psi_a;
			alpha_a*gm1	alpha_a*gm1*M_a*M_a	alpha_a	alpha_a*Psi_a;
				0			0				   1		0	;
				0			0				   0		1];	

	R_a = [	1 	M_a 	0	0;
		   	1  -M_a 	0	0;
			0	0		1	0;
			0	0		0	1];

	s_a = P_a\I_a;

	w_a = R_a*s_a;

	w_p_a_num = w_a(1);
	w_m_a_num = w_a(2);
	w_s_a_num = w_a(3);
	w_z_a_num = w_a(4);

	Res(1) = (w_m_a_num - w_m_a);%(u+c) = prescribed at the inlet
	Res(2) = (w_p_a_num - w_p_a);%(u-c) = prescirbed at the inlet
	Res(3) = (w_s_a_num - w_s_a);%(s) = prescirbed at the inlet
	Res(4) = (w_z_a_num - w_z_a);%(z) = prescribed at the inlet
end%SupersonicBCs()
