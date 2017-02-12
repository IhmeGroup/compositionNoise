function[] = plotInertTransferFunction()
	close all;
	[Nsp, species, a, A, MW] = speciesPropsInert();
	Nmach = 201;
	Npts = 101;
	hh = figure();
	dumb = figure();
	dumber = figure();
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
			[psi_a, ~, ~, ~, ~, gamma] = returnPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
%			[psi_a, ~, ~, ~, ~, ~] = returnPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar);
%			gamma
			GAMMA(i) = gamma;
			Zdp = sqrt(Y*(1-Y));
%			Zdp = 1;
			for j = 1:Nmach
				M_b = (j-1)*0.01;
				Z(i,j) = Y;
				MB(i,j) = M_b;

				T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
				p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
				[psi_b, ~, ~, ~, G(i,j)] = returnPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar,cpfixed);
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

		if (inert == 1)
			figure();
			surface(Z,MB,abs(Zratio),'edgecolor','none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{Ar}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');
			caxis([0,2]);

			figure(dumb)
			plot(Z(:,1), [GAMMA], 'b', 'LineWidth', 2);
			hold on;

			figure();
			surface(Z, MB, PSI, 'edgecolor', 'none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{Ar}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			title('\Psi Ar');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');

			figure();
			surface(Z, MB, G, 'edgecolor', 'none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{Ar}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			title('G Ar');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');
		elseif (inert == 2)
			figure();
			surface(Z,MB,abs(Zratio),'EdgeColor','none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{He}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');
			caxis([0,2]);

			figure(dumb)
			plot(Z(:,1), [GAMMA], 'm', 'LineWidth', 2);
			hold on;

			figure();
			surface(Z, MB, PSI, 'edgecolor', 'none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{Ar}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			title('\Psi He');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');

			figure();
			surface(Z, MB, G, 'edgecolor', 'none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{Ar}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			title('G He');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');
		elseif (inert == 3)
			figure();
			surface(Z,MB,abs(Zratio),'EdgeColor','none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{Kr}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');
			caxis([0,2]);

			figure(dumb)
			plot(Z(:,1), [GAMMA], 'r', 'LineWidth',2);
			hold on;

			figure();
			surface(Z, MB, PSI, 'edgecolor', 'none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{Ar}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			title('\Psi Kr');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');

			figure();
			surface(Z, MB, G, 'edgecolor', 'none');
			shading interp;
			colorbar();
			xlabel('$Y_{\textrm{Ar}}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			title('G Kr');
			set(gca,'FontSize', 14, 'FontName', 'Times');
			colormap('jet');
		end

		figure(hh)
		plot(MB(51,:), [abs(Zratio(51,:))], 'LineWidth',2);
		hold on;

		figure(dumber);
		plot(Z, MYCP./1E4);
		hold on;
		disp('tsup');


	end%inert
	figure(dumber);
	legend('Ar', 'He', 'Kr', 'Ne');
	
	figure(hh);
	plot(MB(51,:), [abs(Pratio(51,:))], 'LineWidth',2);
	plot(MB(51,:), [abs(Sratio(51,:))], 'LineWidth',2);
	xlabel('$M_c$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('Response','FontSize', 14, 'FontName','Times');
	set(gca,'FontSize', 14, 'FontName','Times');
	aa = legend('$\pi^+_c \slash \xi_{a,Ar}$', '$\pi^+_c \slash \xi_{a,He}$', '$\pi^+_c \slash \xi_{a,Kr}$','$\pi^+_c \slash \pi^+_a$', '$\pi^+_c \slash \sigma_a$');
	set(aa, 'Interpreter','LaTeX', 'FontSize', 14, 'FontName','Times', 'Location','NorthWest');

	figure(dumb);
	legend('Argon', 'Helium', 'Krypton');
end%plotInertTransferFunction()
