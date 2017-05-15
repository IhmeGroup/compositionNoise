function[Res] = ShockedBCs(I_a, I_b, w_l, w_r, eta_l, eta_r, flaggo, SPLINES)
%	Boundary conditions for subsonic flows
%
%	Unpack the parameter vector
	[gamma] = returnAmbientState();
%	Pre-compute for convenience	
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	if (flaggo == 0)
		error('flaggo = 0 for shock case');
	elseif (flaggo == 1)
%		Unpack the numerical data from the shock side boundary			
		[w_num_l] = charsFromInvrnts(eta_l, I_a, flaggo, SPLINES);
		w_p_b_p_num = w_num_l(1);
		w_m_b_p_num = w_num_l(2);
		w_s_b_p_num = w_num_l(3);
		w_z_b_p_num = w_num_l(4);

%		Unpack the numerical data from the nozzle exit boundary			
		[w_r_num] = charsFromInvrnts(eta_r, I_b, flaggo, SPLINES);
		w_m_c_num = w_r_num(2);

%		Unpack the supersonic data from the previous solution	
		w_p_b_m	= w_l(1);
		w_m_b_m = w_l(2);
		w_s_b_m = w_l(3);
		w_z_b_m = w_l(4);
		w_z_b_p = w_l(4);

%		Unpack the Mach number from upstream of the shock and then compute the corresponding value downstream from the normal shock relations		
		M_b_m	= SPLINES(2);
		M_b_p	= sqrt((gm1*M_b_m*M_b_m + 2)/(2*gamma*M_b_m*M_b_m - gm1));
		M_c		= SPLINES(3);
		Psibar_b_m	= BaseFlowFromMLVG(M_b_m);
		Psibar_b_p	= BaseFlowFromMLVG(M_b_p);

%		Use the relationships from Luca's JFMR to compute the ideal characteristics on the post shock side
		w_p_b_p = (1 + 2*M_b_p*M_b_p*M_b_m + M_b_m*M_b_m)./(1 + 2*M_b_m*M_b_m*M_b_p + M_b_m*M_b_m)*w_p_b_m ...
				+ (1 - 2*M_b_p*M_b_p*M_b_m + M_b_m*M_b_m)./(1 + 2*M_b_m*M_b_m*M_b_p + M_b_m*M_b_m)*w_m_b_m;
		w_s_b_p	= w_s_b_m + (Psibar_b_p - Psibar_b_m)*w_z_b_m + (gm1*(M_b_m*M_b_m - 1).^2)/(M_b_m*M_b_m*(2 + gm1*M_b_m*M_b_m))*(w_p_b_p + w_m_b_p_num - w_p_b_m - w_m_b_m);

%		Load the residuals to reflect the desired bcs
		Res(1) = w_p_b_p - w_p_b_p_num;
		Res(2) = w_s_b_p - w_s_b_p_num;
		Res(3) = w_z_b_p - w_z_b_p_num;
		Res(4) = 		   w_m_c_num;
	else
		error('Expand flaggo applications');
	end
end
