function[] = plotInertTransferFunction()
	close all;
	addpath('../core');
	addpath('../data');

	fs = 18;
	lw = 4;


	[Nsp, species, a, A, MW] = speciesPropsInert();
	Nmach = 201;
	Npts = 101;
	hh = figure();
	Z = zeros(Npts, Nmach);
	Pratio = zeros(Npts, Nmach);
	Zratio = zeros(Npts, Nmach);
	Sratio = zeros(Npts, Nmach);

	for inert = 1:3
		M_a = 0.0;
		T0 = 300.0;
		p0 = 2.0E6;
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
			[cpfixed]	= returnSpeciesProperties(T0, p0, Yctr, a, A, MW);
			MYCP(i) = cpfixed;
			gamma = cpfixed/(cpfixed - Rbar);
			[psi_a] = returnInertPsi(T0, p0, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
			GAMMA(i) = gamma;
			Zdp = sqrt(Y*(1-Y));
			for j = 1:Nmach
				M_b = (j-1)*0.01;
				Z(i,j) = Y;
				MB(i,j) = M_b;

				T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T0;
				p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p0;
				[psi_b] = returnInertPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar,cpfixed);
				PSI(i,j) = psi_b;
				if (M_b < 1)%subcritical
					Pratio(i,j) = (2*(1 + M_a)*M_b*(2 + (gamma-1.0)*M_b*M_b))/((1+M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
					Sratio(i,j) = ((M_b - M_a)*M_b)/((1+M_b)*(2 + (gamma-1)*M_a*M_b));
					Zratio(i,j) = ((gamma-1)*(psi_b - psi_a)*(2 + (gamma-1)*M_b*M_b)*M_a*M_b)             /((gamma-1)*(1 + M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b)) + ...
	                           (M_b*(2*(psi_a - psi_b) + (gamma-1)*(psi_a*M_b*M_b - psi_b*M_a*M_a)))/((gamma-1)*(1 + M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
					Zratio(i,j) = Zratio(i,j)*Zdp;
				else%supercritical
					Pratio(i,j) = (2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a);
					Sratio(i,j) = (M_b - M_a)/(2*(2 + (gamma-1)*M_a));
					Zratio(i,j) = ((2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a)*psi_a - psi_b)/(2*(gamma-1));
					Zratio(i,j) = Zratio(i,j)*Zdp;
				end%criticality logic
			end%for j = 1:Nmach
		end%for i = 1:Npts

		figure(hh)
		if (inert == 1)
			plot(MB(51,:), [abs(Zratio(51,:))], 'LineWidth',lw);
		elseif (inert == 2)
			plot(MB(51,:), [abs(Zratio(51,:))], '--', 'LineWidth',lw);
		elseif (inert == 3)
			plot(MB(51,:), [abs(Zratio(51,:))], '-.', 'LineWidth',lw);
		end
		hold on;

	end%inert
	figure(hh);
	plot(MB(51,:), [abs(Pratio(51,:))], ':', 'LineWidth',lw);
	plot(MB(51,:), [abs(Sratio(51,:))],  'LineWidth',lw/4);
	plot([1 1], [0, 2], 'k--', 'LineWidth', lw);
	xlabel('$M_c$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
	ylabel('Response','FontSize', fs, 'FontName','Times');
	set(gca,'FontSize', fs, 'FontName','Times');
	aa = legend('$\pi^+_c \slash \xi_{a,Ar}$', '$\pi^+_c \slash \xi_{a,He}$', '$\pi^+_c \slash \xi_{a,Kr}$','$\pi^+_c \slash \pi^+_a$', '$\pi^+_c \slash \sigma_a$');
	set(aa, 'Interpreter','LaTeX', 'FontSize', fs, 'FontName','Times', 'Location','West');

	print -depsc InertTransferFunctionZ50.eps
end%plotInertTransferFunction()
