function[Res] = DMShockedBCs(I_a, I_b, w_l, w_r, eta_l, eta_r, flaggo, SPLINES)
%	Boundary conditions for subsonic flows
%
%	Unpack the parameter vector
	global epsilon;
	global omga;
	Omega = omga;
	[gamma] = returnAmbientState();
%	Pre-compute for convenience	
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	if (flaggo == 0)
		error('flaggo = 0 for shock case');
	elseif (flaggo == 1)
%		Unpack the Mach number from upstream of the shock and then compute the corresponding value downstream from the normal shock relations		
		M_b_m	= SPLINES(2);
		M_b_p	= sqrt((gm1*M_b_m*M_b_m + 2)/(2*gamma*M_b_m*M_b_m - gm1));
		Psibar_b_m	= BaseFlowFromMLVG(M_b_m);
		Psibar_b_p	= BaseFlowFromMLVG(M_b_p);
	elseif (flaggo == 2)
		M_b_m 	=  1.532;
		M_b_p	= 0.716;%sqrt((gm1*M_b_m*M_b_m + 2)/(2*gamma*M_b_m*M_b_m  - gm1))
		[Psibar_b_m	ubar_b_m] = BaseFlowFromMDMSC(eta_l-2*epsilon, M_b_m);
		[Psibar_b_p	ubar_b_p] = BaseFlowFromMDMSC(eta_l, M_b_m);
		Ma1		= MFromEtaDMSC(eta_l - 3*epsilon, SPLINES);
		Ma2 	= MFromEtaDMSC(eta_l - 2*epsilon, SPLINES);
		d_M_b_m_d_xi	 = (Ma2 - Ma1)/epsilon;
		Mb1		= MFromEtaDMSC(eta_l + 1*epsilon, SPLINES);
		Mb2		= MFromEtaDMSC(eta_l + 2*epsilon, SPLINES);
		d_M_b_p_d_xi	= (Mb2 - Mb1)/epsilon;
		c_b_m 	= ubar_b_m./M_b_m;
		c_b_p	= ubar_b_p./M_b_p;
	else
		error('Should have written more code');
	end

	chk = (M_b_p*M_b_p - 1)/(M_b_p*(1 + gm1o2*M_b_p*M_b_p))*d_M_b_p_d_xi - (M_b_m*M_b_m - 1)./(M_b_m*(1 + gm1o2*M_b_m*M_b_m))*d_M_b_m_d_xi;
%	Unpack the numerical data from the shock side boundary			
	[w_num_l] = charsFromInvrnts(eta_l, I_a, flaggo, SPLINES);
	w_p_b_p_num = w_num_l(1);
	w_m_b_p_num = w_num_l(2);
	w_s_b_p_num = w_num_l(3);
	w_z_b_p_num = w_num_l(4);

%	Unpack the numerical data from the nozzle exit boundary			
	[w_r_num] = charsFromInvrnts(eta_r, I_b, flaggo, SPLINES);
	w_m_c_num = w_r_num(2);

%	Unpack the supersonic data from the previous solution	
	w_p_b_m	= w_l(1);
	w_m_b_m = w_l(2);
	w_s_b_m = w_l(3);
	w_z_b_m = w_l(4);
	w_z_b_p = w_l(4);

%	Equations 5.17a-c
	Gamma_p		= (1 - (M_b_p*M_b_p - M_b_m*M_b_m)/(2*M_b_m*M_b_m*M_b_p*M_b_p*(M_b_p*M_b_p - 1))) * (1/M_b_m) *d_M_b_m_d_xi - 2*pi*sqrt(-1)*Omega/(M_b_m*c_b_m);
	Gamma_rho 	= 0.5*(1 + (M_b_m*M_b_m - 1)/(M_b_p*M_b_p - 1))*(1/M_b_m)*d_M_b_m_d_xi - 2*pi*sqrt(-1)*Omega/(M_b_m*c_b_m);
	Gamma_u		= 0.5*(1 + (M_b_m*M_b_m - 1)/(M_b_p*M_b_p - 1))*(1/M_b_m)*d_M_b_m_d_xi - (M_b_m*M_b_m + 1)/2*2*pi*sqrt(-1)*Omega/(M_b_m*c_b_m);

%	Lower half of Table 3
	Lambda 	= (1 - Gamma_rho/Gamma_p)/(M_b_m*(1 + gm1o2*M_b_m*M_b_m));
	delta	= 2*(1 - Gamma_u/Gamma_p)/(1 + gm1o2*M_b_m*M_b_m);
	alpha_p	= Gamma_u/Gamma_p + M_b_m*M_b_p*M_b_p*(1 - delta*(1 - gm1o2*M_b_m));
	alpha_m	= Gamma_u/Gamma_p - M_b_m*M_b_p*M_b_p*(1 - delta*(1 + gm1o2*M_b_m));
	phi		= 0.5*(1 - 1/(M_b_m*M_b_m*M_b_p*M_b_p)*Gamma_rho/Gamma_p);
	Psi_p	= Gamma_u/Gamma_p + M_b_m*M_b_m*M_b_p;
	Psi_m	= Gamma_u/Gamma_p - M_b_m*M_b_m*M_b_p;

%	Upper half of Table 3
	w_p_b_p_w_p_b_m 	= alpha_p/Psi_p;
	w_p_b_p_w_m_b_m		= alpha_m/Psi_p;
	w_p_b_p_w_s_b_m		= M_b_m*M_b_m*M_b_p*M_b_p*delta/Psi_p;
	w_p_b_p_w_m_b_p		= -Psi_m/Psi_p;

	w_s_b_p_w_p_b_m		= phi*(alpha_p/Psi_p - 1) - Lambda*(1 - gm1o2*M_b_m);
	w_s_b_p_w_m_b_m		= phi*(alpha_m/Psi_p - 1) + Lambda*(1 + gm1o2*M_b_m);
	w_s_b_p_w_s_b_m		= 1 + M_b_m*M_b_m*M_b_p*M_b_p*delta/Psi_p + M_b_m*Lambda;
	w_s_b_p_w_m_b_p		= phi*(1 - Psi_m/Psi_p);		

%	Equation 5.18	
	w_p_b_p = w_p_b_p_w_p_b_m*w_p_b_m + w_p_b_p_w_m_b_m*w_m_b_m + w_p_b_p_w_s_b_m*w_s_b_m + w_p_b_p_w_m_b_p*w_m_b_p_num;
%	Equation 5.19	
	w_s_b_p = w_s_b_p_w_p_b_m*w_p_b_m + w_s_b_p_w_m_b_m*w_m_b_m + w_s_b_p_w_s_b_m*w_s_b_m + w_s_b_p_w_m_b_p*w_m_b_p_num;

%	Load the residuals to reflect the desired bcs
	Res(1) = w_p_b_p - w_p_b_p_num;
	Res(2) = w_s_b_p - w_s_b_p_num;
	Res(3) = w_z_b_p - w_z_b_p_num;
	Res(4) = 		   w_m_c_num;
end
