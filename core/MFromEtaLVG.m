function[M] = MFromEtaLVG(eta, SPLINES);
	gamma = returnAmbientState();
	gp1 = gamma + 1;
	gm1 = gamma - 1;
	gp1o2 = gp1/2;
	gm1o2 = gm1/2;
	M_a = SPLINES(1);
	M_b = SPLINES(2);
	M_c = SPLINES(3);
	if (M_c == 0) %subsonic, supersonic, or transonic; all works the same here
		eta_l = sqrt(gp1o2*M_a.^2./(1 + gm1o2*M_a.^2));
		eta_r = sqrt(gp1o2*M_b.^2./(1 + gm1o2*M_b.^2));
		eta = eta*(eta_r - eta_l) + eta_l;
		eta = eta.*eta;
		M	= sqrt(2/gp1*eta/(1 - gm1/gp1*eta));
	else  %shocked
		x_a 	= sqrt(gp1o2*M_a*M_a/(1+gm1o2*M_a*M_a));
		x_bm 	= sqrt(gp1o2*M_b*M_b/(1+gm1o2*M_b*M_b));
		M_p		= sqrt((gm1*M_b*M_b + 2)/(2*gamma*M_b*M_b - gm1));
		x_bp 	= sqrt(gp1o2*M_p*M_p/(1+gm1o2*M_p*M_p));
		x_c		= sqrt(gp1o2*M_c*M_c/(1+gm1o2*M_c*M_c));
		L_1 	= (x_bm - x_a);
		L_2		= (x_bp - x_c);
		L 		= L_1 + L_2;

		y_a		= 0;
		y_b		= L_1/L;
		y_c		= 1;

%		[x_a x_bm x_bp x_c y_a y_b y_c L_1 L_2 L]
%		[M_a M_b M_p M_c gamma]
		if (eta < y_b)
		%	clc
			eta = eta*L+x_a;
			eta = eta.*eta;
			M	= sqrt(2/gp1*eta/(1 - gm1/gp1*eta));
		else
			eta = -(eta - L_1/L)*L + x_bp;
			eta = eta.*eta;
			eta_l = gp1o2*M_a.^2./(1 + gm1o2*M_a.^2);
			eta_r = gp1o2*M_b.^2./(1 + gm1o2*M_b.^2);
			M	= sqrt(2/gp1*eta/(1 - gm1/gp1*eta));

		end
	end
end
