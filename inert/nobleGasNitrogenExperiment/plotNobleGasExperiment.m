function[] = plotInertTransferFunction()
	close all;
	[Nsp, species, a, A, MW] = speciesPropsInert();
	Nsp
	Nmach = 201;
	Npts = 101;
	emdub = [0.1:0.1:0.9]';
%	emdub = [0.001; emdub; 1; 1./flipud(emdub); 1E6]
	emdub = linspace(0.001, .999, 100)';
	emdub = [emdub; 1; 1./flipud(emdub)]
	Nemdub = length(emdub);

	Z 		= zeros(Nemdub, Nmach);
	Pratio 	= zeros(Nemdub, Nmach);
	Zratio 	= zeros(Nemdub, Nmach);
	Sratio 	= zeros(Nemdub, Nmach);
	MWBAR 	= zeros(Nemdub, Nmach);
	MWN2	= 14.0067*2;
	MWO2 	= 15.9994*2;
	global molweight;

	for aa = 1:Nemdub
		gamma = 1.4;
		M_a = 0.0;
		T_a = 300.0;
		p_a = 2.0E6;
		molweight = emdub(aa).*MWN2;
%					Ar		CO		CO2		He		Kr		N2		Ne		O2		H2O	
		MW 		=  [MWN2 MWO2 molweight];
		YN2 	=  [.767 0.233 	0.0000];

		Ynoble 	= [0 0 1];


		for i = 1:1
			Y = 0.5;
			Ylft = (1-(Y-0.01))	*YN2 + (Y-0.01)*Ynoble;
			Yctr = (1-Y)		*YN2 + Y		*Ynoble;
			Yrgt = (1-(Y+0.01))	*YN2 + (Y+0.01)*Ynoble;
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
			GAMMA(aa) = gamma;
%			Zdp = sqrt(Y*(1-Y));
			Zdp = 1;
			for j = 1:Nmach
				M_b = (j-1)*0.01;
				Z(aa,j) = Y;
				MB(aa,j) = M_b;

				T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
				p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
				[psi_b, ~, ~, ~, G(i,j)] = returnPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar,cpfixed);
				PSI(aa,j) = psi_b;
				if (M_b < 1)%subcritical
					Pratio(aa,j) = (2*(1 + M_a)*M_b*(2 + (gamma-1.0)*M_b*M_b))/((1+M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
					Sratio(aa,j) = ((M_b - M_a)*M_b)/((1+M_b)*(2 + (gamma-1)*M_a*M_b));
					Zratio(aa,j) = ((gamma-1)*(psi_b - psi_a)*(2 + (gamma-1)*M_b*M_b)*M_a*M_b)             /((gamma-1)*(1 + M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b)) + ...
	                           (M_b*(2*(psi_a - psi_b) + (gamma-1)*(psi_a*M_b*M_b - psi_b*M_a*M_a)))/((gamma-1)*(1 + M_b)*(M_a + M_b)*(2 + (gamma-1)*M_a*M_b));
					Zratio(aa,j) = Zratio(aa,j)*Zdp;
					MWBAR(aa,j) = molweight;
				else%supercritical
					Pratio(aa,j) = (2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a);
					Sratio(aa,j) = (M_b - M_a)/(2*(2 + (gamma-1)*M_a));
					Zratio(aa,j) = ((2 + (gamma-1)*M_b)/(2 + (gamma-1)*M_a)*psi_a - psi_b)/(2*(gamma-1));
					Zratio(aa,j) = Zratio(aa,j)*Zdp;
					MWBAR(aa,j) = molweight;
				end%criticality logic
			end%for j = 1:Nmach
		end%for i = 1:Npts

		G = G./1E7;
	end%for aa

	MWBAR = MWBAR./MWN2;

	figure();
	surface(MWBAR./(MWBAR + 1), MB, abs(Zratio), 'edgecolor', 'none');
	shading interp;
	hold on;
	plot3([0.5 0.5], [0, 2], [1E6, 1E6], 'w--', 'LineWidth', 4);
	colorbar();
%	xlabel('$MW_x \slash MW_{Ar}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$\beta \slash \left(\beta + 1\right)$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$M_b$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca,'FontSize', 14, 'FontName', 'Times');
	colormap('jet');

	figure();
	surface(MWBAR./(MWBAR + 1), MB, abs(PSI), 'edgecolor', 'none');
	title('Psi');
	shading interp;
	colorbar();
	xlabel('$MW_x \slash MW_{Ar}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$M_b$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca,'FontSize', 14, 'FontName', 'Times');
	colormap('jet');

%	figure();
%	surface(MWBAR./MWN2, MB, abs(PSI), 'edgecolor', 'none');
%	shading interp;
%	colorbar();
%	xlabel('$MW_x \slash MW_{Ar}$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
%	ylabel('$M_b$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
%	set(gca,'FontSize', 14, 'FontName', 'Times');
%	colormap('jet');
end

		
