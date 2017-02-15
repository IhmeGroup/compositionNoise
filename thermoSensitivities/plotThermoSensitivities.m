function[] = plotThermoSensitivities()
	close all;

	lw = 4;
	fs = 22;

	plot_checks = false;%Plot a bunch of things that I can look at to determine whether the grid is sufficiently smooth, etc. 
	NMach = 201;
	addpath('../data');
	addpath('../core');
	fuel = 2;
	Z_st = 0.0551538;
	[Nspecies, species, a, A, MW] = speciesPropsGRI3();
	for run = 1:9
%		First, the chi sweep
		if (run == 1)
			data = load('../data/thermoStudy/CH4_p01_0chi00001tf0300to0300');
		elseif (run == 2)
			data = load('../data/thermoStudy/CH4_p01_0chi00010tf0300to0300');
		elseif (run == 3)
			data = load('../data/thermoStudy/CH4_p01_0chi33.8267tf0300to0300Tst1749');
%		Then, the pre-heat sweep			
		elseif (run == 4)
			data = load('../data/thermoStudy/CH4_p01_0chi00001tf0400to0400');
		elseif (run == 5)
			data = load('../data/thermoStudy/CH4_p01_0chi00001tf0600to0600');
		elseif (run == 6)
			data = load('../data/thermoStudy/CH4_p01_0chi00001tf0900to0900');
%		Lastly, the pressure sweep			
		elseif (run == 7)
			data = load('../data/thermoStudy/CH4_p16_0chi00001tf0300to0300');
		elseif (run == 8)
			data = load('../data/thermoStudy/CH4_p40_0chi00001tf0300to0300');
		elseif (run == 9)
			data = load('../data/thermoStudy/CH4_p64_0chi00001tf0300to0300');
		end
		[~, Npts] = size(data);
		hh = figure();
		M_a = 0.0;
		p_a = 1.0E6;
		MB = zeros(NMach,Npts)';
		ZEE = zeros(NMach,Npts)';
		Pratio = zeros(NMach,Npts)';
		Sratio = zeros(NMach,Npts)';
		Zratio = zeros(NMach,Npts)';
		TA = zeros(Npts,1);
		PSIA = zeros(Npts,1);
		ZA = zeros(Npts,1);
		CP = zeros(Npts,1);
		HA = zeros(Npts,1);
		SA = zeros(Npts, 1);
		GA = zeros(Npts,1);
		YA = zeros(Npts,Nspecies);
		for i = 1:Npts
			if (run == 4)
				T_a = 2100.0;
			else
				T_a  = data(2,i);
			end
			TA(i) = T_a;
			if (i == 1)
				Ylft = data(3:Nspecies+2,i+1);
				Yctr = data(3:Nspecies+2,i);
				Yrgt = data(3:Nspecies+2,i+2);
				Zlft = data(1,i+1);
				Zctr = data(1,i);
				Zrgt = data(1,i+2);
			elseif (i == Npts)
				Ylft = data(3:Nspecies+2,i-1);
				Yctr = data(3:Nspecies+2,i);
				Yrgt = data(3:Nspecies+2,i-2);
				Zlft = data(1,i-1);
				Zctr = data(1,i);
				Zrgt = data(1,i-2);
			else
				Ylft = data(3:Nspecies+2,i-1);
				Yctr = data(3:Nspecies+2,i);
				Yrgt = data(3:Nspecies+2,i+1);
				Zlft = data(1,i-1);
				Zctr = data(1,i);
				Zrgt = data(1,i+1);
			end
			MWbar = 0;
			for ll = 1:Nspecies
				MWbar = MWbar + Yctr(ll)/MW(ll);
			end
			MWbar = 1/MWbar;
			Rbar = 8.31446E7/MWbar;
			ZA(i) = Zctr./(Zctr + Z_st);
			YA(i,:) = Yctr;
			[cpfixed] = returnSpeciesProperties(T_a, p_a, Yctr, a, A, MW);
			[psi_a, CP(i), HA(i), SA(i), GA(i), gamma] = returnInertPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
			PSIA(i) = psi_a;
			for j = 1:NMach
				M_b = (j-1)*0.01;
				T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
				p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
				MB(i,j) = M_b;
				phi = Zctr./(1-Zctr).*(1-Z_st)./Z_st;
%				ZEE(i,j) = Zctr./(Zctr + Z_st);
				ZEE(i,j) = phi/(1+phi);
				psi_b = returnInertPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
				if (M_b < 1)%subcritical
					Pratio(i,j) = (2*(1 + M_a)*M_b*(2 + (gamma-1.0)*M_b*M_b))/((1+M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
					Sratio(i,j) = ((M_b - M_a)*M_b)/((1+M_b)*(2 + (gamma-1)*M_a*M_b));
					Zratio(i,j) = ((gamma-1)*(psi_b - psi_a)*(2 + (gamma-1)*M_b*M_b)*M_a*M_b)/((gamma-1)*(1+M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b)) + ...
		                           (M_b*(2*(psi_a - psi_b) + (gamma-1)*(psi_a*M_b*M_b - psi_b*M_a*M_a)))/((gamma-1)*(1 + M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
				else%supercritical
					Pratio(i,j) = (2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a);
					Sratio(i,j) = (M_b - M_a)/(2*(2 + (gamma-1)*M_a));
					Zratio(i,j) = ((2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a)*psi_a - psi_b)/(2*(gamma-1));
				end%criticality logic
			end%for j = 1:M_b
	
		end%i = 1:Npts

%		Interpolate the ends b/c the data is garbage
%		Zratio(1,:) = (Zratio(3,:) - Zratio(2,:))/(ZA(3) - ZA(2)) * (ZA(1) - ZA(2)) + Zratio(2,:);
%		Zratio(Npts,:) = (Zratio(Npts-2,:) - Zratio(Npts-1,:))/(ZA(Npts-2) - ZA(Npts-1)) * (ZA(Npts-1) - ZA(Npts)) + Zratio(Npts-1,:);

		figure(hh)
		surface(ZEE, MB, log10(abs(Zratio./Sratio)), 'EdgeColor','None');
		shading('interp');
	%	surface(log10(abs(Zratio)), 'EdgeColor','None');
		xlabel('$Z^*$ [-]', 'Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
		ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
		set(gca,'FontSize', fs, 'FontName','Times');
		gg = colorbar('northoutside');
		set(gg, 'FontSize', fs, 'FontName', 'Times');
		colormap('jet');
		caxis([-1,2.5])
%		xlim([0, 1/(1 + Z_st)]);
		xlim([0, 1]);
		hold on;
		plot3([0 1],[1 1], [1E9, 1E9], 'w--', 'LineWidth', lw); 
		if (run == 1)
			print -djpeg CH4_chi_01.jpg
		elseif (run == 2)
			print -djpeg CH4_chi_10.jpg
		elseif (run == 3)
			print -djpeg CH4_chi_33.jpg
		elseif (run == 4)
			print -djpeg CH4_T_400.jpg
		elseif (run == 5)
			print -djpeg CH4_T_600.jpg
		elseif (run == 6)
			print -djpeg CH4_T_900.jpg
		elseif (run == 7)
			print -djpeg CH4_p_16.jpg
		elseif (run == 8)
			print -djpeg CH4_p_40.jpg
		elseif (run == 9)
			print -djpeg CH4_p_64.jpg
		end

		if (plot_checks)
			figure();
			plot(TA, PSIA, '-o');
			xlabel('$T_a$ [K]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$\psi_a$ [-]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	
			figure();
			plot(ZA, PSIA, '-o');
			xlabel('$Z$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$\psi_a$ [-]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

			figure();
			plot(ZA, GA, '-o');
			xlabel('$Z$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$g_a$ [-]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

			figure();
			plot(ZA, SA, '-o');
			xlabel('$Z$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$s_a$ [-]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

			figure();
			plot(ZA, HA, '-o');
			xlabel('$Z$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$h_a$ [-]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

			figure();
			plot(ZA, CP, '-o');
			xlabel('$Z$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$c_p$ [-]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

			figure();
			plot(ZA, TA, '-o');
			xlabel('$Z$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$T_a$ [-]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

			deltaY = 5;
			zz = 0;
%			for zz = 1:deltaY:Nspecies
			while (zz < Nspecies/deltaY)
				zz = zz + 1;
				figure();
				size(YA(:,(zz-1)*deltaY+1:min(zz*deltaY,Nspecies)))
				semilogy(ZA, YA(:,(zz-1)*deltaY+1:min(zz*deltaY,Nspecies)), '-o');
				xlabel('$Z$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
				ylabel('$\Sigma Y_a$ [-]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
				legend(species((zz-1)*deltaY+1:min(zz*deltaY,Nspecies)));
			end
		end%(plot_checks)
end%plotFuelTransferFunctions()
