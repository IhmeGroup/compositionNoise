function[] = plotPressureSensitivity()

	fs = 22;



	[Nsp, species, a, A, MW] = speciesPropsInert();
%	gamma = 1.4;
	PA 		= zeros(101,201);
	MB 		= zeros(101,201);
	Zratio 	= zeros(101,201);
	M_a = 0.0;

	MW 		=  [39.948 28.0101 44.0095  4.0026  83.798  28.0134	20.1791 31.9988 18.0153]; 
	Yair 	=  [0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.7670 	0.0000 	0.2330	0.0000];
	Y = 0.5;
	Yfuel= [0.0000 	0.0000 	0.0000 	1.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000];%Helium
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
	for i = 1:101
%		T_a = 290 + i*10;
		T_a = 300.0;
		p_a = 9.0E5 + i*1.0E5;
%				Ar		CO		CO2		He		Kr		N2		Ne		O2		H2O	
%		Ylft = [0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000	0.0000];
%		Yrgt = [0.0000 	0.0000 	0.0000 	1.0000 	0.0000 	0.0000 	0.0000 	0.0000 	0.0000];
%		Yctr = 0.5*(Ylft + Yrgt);
%		Zlft = 0.0;
%		Zctr = 0.5;
%		Zrgt = 1.0;

%		gamma = 1.4;
%		psi_a = 					 returnInertPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar);
		[cpfixed]	= returnSpeciesProperties(T_a, p_a, Yctr, a, A, MW);
		[psi_a, ~, ~, ~, ~, gamma] = returnInertPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
		for j = 1:201
			M_b = (j-1)*0.01;
			T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
			p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
			MB(i,j) = M_b;
			PA(i,j) = p_a;
			psi_b = returnInertPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar,cpfixed);
			if (M_b < 1)%subcritical
				Zratio(i,j) = ((gamma-1)*(psi_b - psi_a)*(2 + (gamma-1)*M_b*M_b)*M_a*M_b)/((gamma-1)*(1+M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b)) + ...
	                            (M_b*(2*(psi_a - psi_b) + (gamma-1)*(psi_a*M_b*M_b - psi_b*M_a*M_a)))/((gamma-1)*(1 + M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
			else%supercritical
				Zratio(i,j) = ((2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a)*psi_a - psi_b)/(2*(gamma-1));
			end%criticality logic
		end%for M_b
	end%for T_a
	Zratio = abs(Zratio);

	hh = figure();
	surface(MB, PA./1.0E6, Zratio,'edgecolor','none');
	shading interp;
	colormap('jet');
	h = colorbar('northoutside');
%	set(h, 'FontSize', fs, 'FontName' 'Times');
	xlabel('$M_c$','Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
	ylabel('$p_a$ [bar]', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName','Times');
	set(gca,'FontSize', fs, 'FontName','Times');
	xlim([0, 2]);
	ylim([1.0E6, 1.0E7]./1.0E6);

	print -djpeg InertPressureSensitivity.jpg

end
