function[plus_d, minus_d, sigma_d, xi_d] = compactNoise(type, M_a, M_b, M_d)
%	Here, 	a refers to the nozzle inlet condition, 
%			b refers to the supersonic state immediately before the shock
%			c refers to the subsonic state immediately after the shock
% 	and 	d refers to the nozzle exit condition
	disp('hello');

	plus_a 	= 0;
	minus_d = 0;
	sigma_a = 0;
	xi_a 	= 0;
	if (type == 1)		plus_a 	= 1; 	%Downstream acoustic forcing
	elseif (type == 2) 	minus_d	= 1; 		%Upstream acoustic forcing
	elseif (type == 3) 	sigma_a = 1;	%Entropy forcing
	elseif (type == 4) 	xi_a 	= 1;	%Compositional forcing
	end


	addpath('../core');
	addpath('../speciesProps');
	global Nspecies species a A MW Hover;
	[Nspecies, species, a, A, MW, Hover] = speciesPropsCH4();
	global mechanism;
	mechanism = 2;
	loadSpeciesData();
	
%	Prescribe the stagnation condition	
	T0 	= 1600;% K
	p0	= 4.5E5;% Pa
	Zbar = 0.03;

	gamma = returnGamma(T0, p0, Zbar);
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;

%	Apply the normal shock relations to determine M_c
	M_c = sqrt((gm1*M_b*M_b + 2)/(2*gamma*M_b - gm1));


	disp('ello');
%	Compute some base state stuff	
	T_a		= (1 + gm1o2*M_a*M_a).^(-1)*T0;
	p_a		= (1 + gm1o2*M_b*M_b).^(-gamma/gm1)*p0;
	Psi_a = returnPsi(T_a, p_a, Zbar);

	T_b		= (1 + gm1o2*M_b*M_b).^(-1)*T0;
	p_b		= (1 + gm1o2*M_b*M_b).^(-gamma/gm1)*p0;
	Psi_b = returnPsi(T_b, p_b, Zbar);

%	Adjust the stagnation pressure to the post shock state
	p0 = (gp1*M_b*M_b/(gm1*M_b*M_b + 2)).^(gamma/gm1)*(gp1/(2*gamma*M_b*M_b - gm1)).^(1/gm1);

	T_c		= (1 + gm1o2*M_c*M_c).^(-1)*T0;
	p_c		= (1 + gm1o2*M_c*M_c).^(-gamma/gm1)*p0;
	Psi_c = returnPsi(T_c, p_c, Zbar);

	T_d		= (1 + gm1o2*M_d*M_d).^(-1)*T0;
	p_d		= (1 + gm1o2*M_d*M_d).^(-gamma/gm1)*p0;
	Psi_d = returnPsi(T_d, p_d, Zbar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Part 1 - Nozzle inlet to shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Hard code the transfer functions
	plus_b_over_plus_a 		= (2 + gm1*M_b)/(2+gm1*M_a);
	plus_b_over_sigma_a		= (M_b - M_a)/(2*(2 + gm1*M_a));
	plus_b_over_xi_a		= 1/(2*gm1)*(-Psi_b + (2 + gm1*M_b)/(2 + gm1*M_a)*Psi_a);

	minus_b_over_plus_a		= (2 - gm1*M_b)/(2 + gm1*M_a);
	minus_b_over_sigma_a	= -(M_a + M_b)/(2*(2 + gm1*M_a));
	minus_b_over_xi_a		= -1/(2*gm1)*(Psi_b + (gm1*M_b - 2)/(gm1*M_a + 2)*Psi_a);

%	Apply the transfer functions to the disturbances	
	plus_b 	= plus_b_over_plus_a*plus_a 	+ plus_b_over_sigma_a*sigma_a 	+ plus_b_over_xi_a*xi_a;
	sigma_b	= sigma_a;
	minus_b	= minus_b_over_plus_a*plus_a 	+ minus_b_over_sigma_a*sigma_a 	+ minus_b_over_xi_a*xi_a;
	xi_b	= xi_a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Part 2 - Across the shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Hard code the transfer functions
	plus_c_over_plus_b		= (1 + 2*M_c*M_c*M_b + M_b*M_b)/(1 + 2*M_b*M_b*M_c + M_b*M_b);
	plus_c_over_minus_b 	= (1 - 2*M_c*M_c*M_b + M_b*M_b)/(1 + 2*M_b*M_b*M_c + M_b*M_b);
	sigma_c_over_xi_b		= -Psi_b;
	sigma_c_over_xi_c		= Psi_c;
	sigma_c_over_plus_c		= gm1*(M_b*M_b - 1)*(M_b*M_b - 1)/(M_b*M_b*(2 + gm1*M_b*M_b));
	sigma_c_over_minus_c	= gm1*(M_b*M_b - 1)*(M_b*M_b - 1)/(M_b*M_b*(2 + gm1*M_b*M_b));
	sigma_c_over_plus_b		=-gm1*(M_b*M_b - 1)*(M_b*M_b - 1)/(M_b*M_b*(2 + gm1*M_b*M_b));
	sigma_c_over_minus_b	=-gm1*(M_b*M_b - 1)*(M_b*M_b - 1)/(M_b*M_b*(2 + gm1*M_b*M_b));

%	Hard code the transfer functions
	plus_d_over_plus_c 		= 2*(1 + M_c)*M_b*(2 + gm1*M_d*M_d)/((1 + M_d)*(M_c + M_d)*(2 + gm1*M_c*M_d));
	plus_d_over_sigma_c		= (M_d - M_c)*M_d/((1+M_d)*(2 + gm1*M_c*M_d));
	plus_d_over_xi_c		= gm1*(Psi_d - Psi_c)*(2 + gm1*M_d*M_d)*M_c*M_d/(gm1*(1 + M_d)*(M_c+M_d)*(2 + gm1*M_c*M_d)) ...
							+ M_d*(2*(Psi_c - Psi_d) + gm1*(Psi_c*M_d*M_d - Psi_d*M_c*M_c))/(gm1*(1 + M_d)*(M_c + M_d)*(2 + gm1*M_c*M_d));
	minus_c_over_plus_c		= (M_d - M_c)*(1 - M_c)*(1 - gm1o2*M_c*M_d)/((1 - M_c)*(M_c + M_d)*(1 + gm1o2*M_c*M_d));
	minus_c_over_sigma_c	= -(M_d - M_c)*M_c/((1 - M_c)*(2 + gm1*M_c*M_d));
	minus_c_over_xi_c		= gm1*(Psi_c - Psi_d)*(2 + gm1*M_c*M_c)*M_c*M_d/(gm1*(M_c - 1)*(M_c + M_d)*(2 + gm1*M_c*M_d)) ...
							+ (M_c*(2*(Psi_c - Psi_d) + gm1*(Psi_c*M_d*M_d - Psi_d*M_c*M_c)))/(gm1*(M_c - 1)*(M_c + M_d)*(2 + gm1*M_c*M_d));

%	Apply the transfer functions to the disturbances
	plus_c 	= plus_c_over_plus_b*plus_b + plus_c_over_minus_b*minus_b;
	xi_c	= xi_b;

%	De-coupling the subsonic charcteristic from entropy generation at the shock
	a 		= sigma_b + sigma_c_over_xi_b*xi_b + sigma_c_over_xi_c*xi_c + sigma_c_over_plus_c*plus_c + sigma_c_over_plus_b*plus_b + sigma_c_over_minus_b*minus_b;
	b		= sigma_c_over_minus_c;
	alpha	= minus_c_over_plus_c*plus_c + minus_c_over_xi_c*xi_c;
	beta	= minus_c_over_sigma_c;

%	Keep applying the transfer functions to the disturbances	
	sigma_c	= (a + alpha*b)/(1 - beta*b);
	minus_c	= (alpha + a*b)/(1 - beta*b);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Part 3 - Post-Shock to Exit (Subsonic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Apply the transfer functions to the disturbances
	disp('ho');
	plus_d	= plus_d_over_plus_c*plus_c + plus_d_over_sigma_c*sigma_c + plus_d_over_xi_c*xi_c;
	sigma_d	= sigma_c;
	xi_d	= xi_c;
	minus_d	= 0; %Prescribed
end
