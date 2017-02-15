function[] = plotInertShockTransferFunctions()
	close all;

	fs = 22;
	lw = 4;


	addpath('../data');
	addpath('../core');

	[Nsp, species, a, A, MW] = speciesPropsInert();
	Nmach = 301;
	Npts = 101;
	Z = zeros(Npts, Nmach);
	Pratio = zeros(Npts, Nmach);
	Zratio = zeros(Npts, Nmach);
	Sratio = zeros(Npts, Nmach);

	for inert = 1:3
		gamma = 1.4;
		M_a = 0.0;
		T_a = 300.0;
		p_a = 2.0E6;
%					Ar		CO		CO2		He		Kr		N2		Ne		O2		H2O	
		MW 		=  [39.948 28.0101 44.0095  4.0026  83.798  28.0134	20.1791 31.9988 18.0153]; 
		Yair 	=  [0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.7670 	0.0000 	0.2330	0.0000];

		if (inert == 1)%Argon
			Yfuel= [1.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000];
		elseif (inert == 2) %Helium
			Yfuel= [0.0000 	0.0000 	0.0000 	1.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000];
		elseif (inert == 3) %Krypton
			Yfuel= [0.0000 	0.0000 	0.0000 	0.0000 	1.0000 	0.0000 	0.0000 	0.0000 	0.0000];
		elseif (inert == 4) %Neon
			Yfuel= [0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	1.0000 	0.0000 	0.0000];
		end


		for i = 1:Npts
			Y = (i-1)*0.01;
			Ylft = (1-(Y-0.01))	*Yair + (Y-0.01)*Yfuel;
			Yctr = (1-Y)		*Yair + Y		*Yfuel;
			Yrgt = (1-(Y+0.01))	*Yair + (Y+0.01)*Yfuel;
			Zlft = (Y-0.01);
			Zctr = Y;
			Zrgt = (Y+0.01);
			MWbar = 0;
			for ll = 1:Nsp
				MWbar = MWbar + Yctr(ll)/MW(ll);
			end
			MWbar = 1/MWbar;
			Rbar = 8.31446E7/MWbar;
			MYMW(i) = MWbar;
			MYR(i) = Rbar;
			[cpfixed]	= returnSpeciesProperties(T_a, p_a, Yctr, a, A, MW);
			MYCP(i) = cpfixed;
			[psi_a, ~, ~, ~, ~, gamma] = returnInertPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
%			[psi_a, ~, ~, ~, ~, ~] = returnInertPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar);
%			gamma
			GAMMA(i) = gamma;
			gm1 = gamma - 1;
			gp1 = gamma + 1;
			for j = 1:Nmach
				M_b = 1 + (j-1)*0.01;
				Z(i,j) = Y;
				MB(i,j) = M_b;

				T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
				p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
				[psi_b, ~, ~, ~, G(i,j)] = returnInertPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar,cpfixed);

			
				M_c = (gm1*M_b^2 + 2)/(2*gamma*M_b^2 - gm1);
				T_c = T_b*(2*gamma*M_b^2 - gm1)*(gm1*M_b^2 + 2)/(gp1^2*M_b^2);
				p_c = gp1*M_b^2/(gm1*M_b^2 + 2);
%				MC(i,j) = M_C
				[psi_c, ~, ~, ~, ~] = returnInertPsi(T_c, p_c, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);

				PSI(i,j) = psi_b;
%				if (M_b < 1)%subcritical
%					Pratio(i,j) = (2*(1 + M_a)*M_b*(2 + (gamma-1.0)*M_b*M_b))/((1+M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
%					Sratio(i,j) = ((M_b - M_a)*M_b)/((1+M_b)*(2 + (gamma-1)*M_a*M_b));
%					Zratio(i,j) = ((gamma-1)*(psi_b - psi_a)*(2 + (gamma-1)*M_b*M_b)*M_a*M_b)             /((gamma-1)*(1 + M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b)) + ...
%	                           (M_b*(2*(psi_a - psi_b) + (gamma-1)*(psi_a*M_b*M_b - psi_b*M_a*M_a)))/((gamma-1)*(1 + M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
%				else%supercritical
%				Pratio(i,j) = (2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a);
%				Sratio(i,j) = (M_b - M_a)/(2*(2 + (gamma-1)*M_a));
%				Zratio(i,j) = ((2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a)*psi_a - psi_b)/(2*(gamma-1));
%				end%criticality logic

				Pbminus = -1/(2*gm1)*(psi_b + (gm1*M_b - 2)/(gm1*M_a + 2)*psi_a);
				Pbplus 	= 1/(2*gm1)*(-psi_b + (2 + gm1*M_b)/(2+gm1*M_a)*psi_a);
	
				Pcplus = (1 + 2*M_c^2*M_b + M_b^2)/(1+2*M_b^2*M_c + M_b^2)*Pbplus + (1-2*M_c^2*M_b + M_b^2)/(1 + 2*M_b^2*M_c + M_b^2)*Pbminus;

				Zratio(i,j) = Pcplus;

			end%for j = 1:Nmach
		end%for i = 1:Npts

		figure();
		surface(Z,MB,log10(abs(Zratio)),'edgecolor','none');
		shading interp;
		gg = colorbar('northoutside');
		set(gg, 'FontSize', fs, 'FontName', 'Times');
		ylabel('$M_{b-\epsilon}$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
		set(gca,'FontSize', fs, 'FontName', 'Times');
		colormap('jet');
		caxis([-1,1]);
		if (inert == 1)
			xlabel('$Y_{\textrm{Ar}}$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
			print -djpeg ShockResultsAr.jpg
		elseif (inert == 2)
			xlabel('$Y_{\textrm{He}}$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
			print -djpeg ShockResultsHe.jpg
		elseif (inert == 3)
			xlabel('$Y_{\textrm{Kr}}$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
			print -djpeg ShockResultsKr.jpg
		end
	end%inert
end%plotInertTransferFunction()
