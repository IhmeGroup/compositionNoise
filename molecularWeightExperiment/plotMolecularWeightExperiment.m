function[] = plotMolecularWeightExperiment()
	close all;

	fs = 22;
	lw = 4;

	addpath('../core');
	addpath('../data');
	[Nsp, species, a, A, MW] = speciesPropsMW();
	Nmach = 201;
	Npts = 101;

	M_a = 0.0;
	T_a = 300.0;
	p_a = 2.0E6;
	MW = 		[2*14.00670 39.94800 40.06500];
	YN2 	= [1 0 0];

	for pass = 1:2
		if (pass == 1)
			Ytest 	=  [1.0000 	1.0000 0.0000];
		elseif (pass == 2)
			Ytest 	=  [1.0000 	0.0000 1.0000];
		else
			error(0);
		end

		Z 		= zeros(Npts, Nmach);
		Pratio 	= zeros(Npts, Nmach);
		Zratio 	= zeros(Npts, Nmach);
		Sratio 	= zeros(Npts, Nmach);

		for i = 1:Npts
			Y = (i-1)*0.01;
			Ylft = (1-(Y-0.01))	*Ytest + (Y-0.01)*YN2;
			Yctr = (1-Y)		*Ytest + Y		*YN2;
			Yrgt = (1-(Y+0.01))	*Ytest + (Y+0.01)*YN2;
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
			GAMMA(i) = gamma;
%			Zdp = sqrt(Y*(1-Y));
			Zdp = 1;
			for j = 1:Nmach
				M_b = (j-1)*0.01;
				Z(i,j) = Y;
				MB(i,j) = M_b;

				T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
				p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
				[psi_b, ~, ~, ~, G(i,j)] = returnInertPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar,cpfixed);
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
		G = G./1E7;
	
		Zratio(end,:) = (Zratio(end-1,:) - Zratio(end-2,:))/(Z(end-1,:) - Z(end-2,:))*(Z(end,:) - Z(end-1,:)) + Zratio(end-1,:);
		Zratio(1,:) = (Zratio(3,:) - Zratio(2,:))/(Z(3,:) - Z(2,:))*(Z(1,:) - Z(2,:)) + Zratio(2,:);


		figure();
		surface(Z, MB, abs(Zratio), 'edgecolor', 'none');
		shading interp;
		hold on;
		plot3([0. 1.], [1, 1], [1E6, 1E6], 'w--', 'LineWidth', lw);
		gg = colorbar('northoutside');
		set(gg, 'FontSize', fs, 'FontName', 'Times');
%		xlabel('$Z \slash \left(Z + Z_{st}\right)$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
		xlabel('$Z$ [-]', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
		ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
		set(gca,'FontSize', fs, 'FontName', 'Times');
		colormap('jet');

		if (pass == 1)
			print -djpeg ArgonSweep.jpg
		elseif (pass == 2)
			print -djpeg C3H4Sweep.jpg
		end

	end%for pass
end%function
