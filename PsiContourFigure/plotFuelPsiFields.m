function[] = plotFuelPsiFields(fuel)
	close all;
	NMach = 201;

	addpath('../data');
	addpath('../core');
	if (fuel == 1)
		Z_st = 0.0627964;
		[Nspecies, species, a, A, MW] = speciesPropsC12H26();
	elseif (fuel == 2)
		Z_st = 0.0551538;
		[Nspecies, species, a, A, MW] = speciesPropsCH4();
	elseif (fuel == 3)
		Z_st = 0.0285207;
		[Nspecies, species, a, A, MW] = speciesPropsH2();
	else
		error('Wrong fuel type selected');
	end
	data = loadFuelData(fuel);
	[~, Npts] = size(data);


	hh = figure();
	M_a = 0.0;
	p_a = 1.0E6;
	MB = zeros(NMach,Npts)';
	ZEE = zeros(NMach,Npts)';
	ZA = zeros(Npts,1);
	PSIB = zeros(Npts,1);
	for i = 1:Npts
		T_a  = data(2,i);
		TA(i) = T_a;
		if (i == 1)
			Ylft = data(3:Nspecies+2,i);
			Yctr = data(3:Nspecies+2,i);
			Yrgt = data(3:Nspecies+2,i+1);
			Zlft = data(1,i);
			Zctr = data(1,i);
			Zrgt = data(1,i+1);
		elseif (i == Npts)
			Ylft = data(3:Nspecies+2,i-1);
			Yctr = data(3:Nspecies+2,i);
			Yrgt = data(3:Nspecies+2,i);
			Zlft = data(1,i-1);
			Zctr = data(1,i);
			Zrgt = data(1,i);
		else
			Ylft = data(3:Nspecies+2,i-1);
			Yctr = data(3:Nspecies+2,i);
			Yrgt = data(3:Nspecies+2,i+1);
			Zlft = data(1,i-1);
			Zctr = data(1,i);
			Zrgt = data(1,i+1);
		end
%		Account for variable gamma
		MWbar = dot(MW, Yctr);
		Rbar = 8.31446E7/MWbar;
		ZA(i) = Zctr;
		[cpfixed, ~, ~, ~] = returnSpeciesProperties(T_a, p_a, Yctr, a, A, MW);
%		[psi_a,	gamma] = returnPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
		[psi_a] = returnPsi(T_a, p_a, Zctr);
		gamma = cpfixed/(cpfixed - Rbar);
		for j = 1:NMach
			M_b = (j-1)*0.01;
			T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
			p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
			MB(i,j) = M_b;
			ZEE(i,j) = Zctr./(Zctr + Z_st);
%			psi_b = returnPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
			psi_b = returnPsi(T_b, p_b, Zctr);
			PSIB(i,j) = psi_b;
		end%for j = 1:M_b
	end%i = 1:Npts

%	Interpolate the ends b/c the data is garbage
	PSIB(1,:) = (PSIB(3,:) - PSIB(2,:))/(ZA(3) - ZA(2)) * (ZA(1) - ZA(2)) + PSIB(2,:);
	PSIB(Npts,:) = (PSIB(Npts-2,:) - PSIB(Npts-1,:))/(ZA(Npts-2) - ZA(Npts-1)) * (ZA(Npts-1) - ZA(Npts)) + PSIB(Npts-1,:);

	figure(hh)
	surface(ZEE, MB, PSIB, 'EdgeColor','None');
	shading('interp');
	xlabel('$Z \slash (Z + Z_{st})$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$M_b$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	title('$\psi$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca,'FontSize', 14, 'FontName','Times');
	colorbar();
	colormap('jet');
	xlim([0, 1/(1 + Z_st)]);

	h2 = figure();
	dPSIBdMach = zeros(Npts, NMach);
	for i = 1:Npts
		dPSIBdMach(i,1) = (PSIB(i,2) - PSIB(i,1))/(MB(i,2) - MB(i,1));
		for j = 2:NMach - 1
			dPSIBdMach(i,j) = (PSIB(i,j+1) - PSIB(i,j-1))/(MB(i,j+1) - MB(i,j-1));
		end%for j
		dPSIBdMach(i,NMach) = (PSIB(i,NMach) - PSIB(i,NMach-1))/(MB(i,NMach) - MB(i,NMach-1));
	end
	surface(ZEE, MB, dPSIBdMach, 'EdgeColor', 'None');
	shading('interp');
	xlabel('$Z \slash (Z + Z_{st})$ [-]', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$M_b$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	title('$\partial \psi \slash \partial M_b$','Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca,'FontSize', 14, 'FontName','Times');
	colorbar();
	colormap('jet');
	xlim([0, 1/(1 + Z_st)]);
end%plotFuelPsiFields
