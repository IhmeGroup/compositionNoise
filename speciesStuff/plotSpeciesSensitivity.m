function[] = plotSpeciesSensitivity()
	close all;

	fs = 22;
	lw = 4;


%	Flags
	mixture_stuff = false
	bar_graph = true;

%	Parameters
	NMach = 201;
	T_a = 300;
	p_a = 1.0E6;

%	Load flamelet data
	data = load('../data/lowStrain/lowStrain.CH4');
	Z_st = 0.0551538;
	[~, Npts] = size(data);
	[Nspecies, species, a, A, MW] = speciesPropsCH4();

%	allocate
	PSI = zeros(NMach,1);
	MB = zeros(NMach,1);
	G = zeros(NMach,1);
	GI = zeros(NMach, Nspecies);
	CP = zeros(NMach, 1);
	CPI = zeros(NMach, Nspecies);
	PSII  = zeros(NMach, Nspecies);

%	indices = [14, 26, 81];
	indices= [9, 26, 99];
	for kase = 1:3
		kase
		condlft = data(:, indices(kase) - 1);
		condctr = data(:, indices(kase));
		condrgt = data(:, indices(kase) + 1);
		T_a = condctr(2);
		cond = data(:,indices(kase));
		Ylft = condlft(3:end);
		Yctr = condctr(3:end);
		Yrgt = condrgt(3:end);
		Zlft = condlft(1);
		Zctr = condctr(1);
		Zrgt = condrgt(1);
		MWbar = 0;
		for ll = 1:Nspecies
			MWbar = MWbar + Yctr(ll)/MW(ll);
		end
		MWbar = 1/MWbar;
		Rbar = 8.31446E7/MWbar;
		[cpfixed]	= returnSpeciesProperties(T_a, p_a, Yctr, a, A, MW);
		[psi_a, ~, ~, ~, ~, gamma] = returnPsi(T_a, p_a, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
		for i = 1:NMach
			M_b = (i-1)*0.01;
			MB(i) = M_b;
			T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
			p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
			[PSI(i), CP(i), ~, ~, G(i)] = returnPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed);
			for j = 1:Nspecies
				Y = zeros(Nspecies,1);
				Y(j) = 1.0;
				[CPI(i,j), ~, ~, GI(i,j)] = returnSpeciesProperties(T_b, p_b, Y, a, A, MW);
%				[PSII(i,:)] = returnSpeciesPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW);
			end%j = 1:Nspecies
		end%i = 1:NMach	
%		1 erg/gram is 1E-7 kJ/kg
		CP = CP*1E-7;
		CPI = CPI*1E-7;
		G = G*1E-10;
		GI = GI*1E-10;

		if (bar_graph)
			ll = figure();
%			set(ll, 'Position', [0 0 650 450]);
			for iii = 1:1:3
				M_b = (iii-1);
				T_b = (1 + (gamma-1)/2*M_b*M_b)^(-1)*T_a;
				p_b = (1 + (gamma-1)/2*M_b*M_b)^(-gamma/(gamma-1))*p_a;
				disp('oats');
				psii = returnSpeciesPsi(T_b, p_b, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW);
%				if (iii == 1)
%					[~, I] = sortrows(abs(psii));
%					I = flipud(I);
%				end

				I = [16, 03, 30, 29, 17, 05, 04, 10];
				Nspecial = length(I');
				bardata = zeros(Nspecial,1);
				barlabels = {};
				for i = 1:Nspecial 
					bardata(i) = psii(I(i))./1E3;
					barlabels(i) = species(I(i));
				end
				if (iii == 1)
%					aaa = bar(log10(abs(bardata)));
					aaa = bar(abs(bardata));
					XX = get(aaa, 'XData');
					set(aaa, 'FaceColor', 'k', 'EdgeColor', 'k', 'XData', XX + -1/4, 'FaceAlpha', 1, 'LineWidth', 2, 'BarWidth', 0.25, 'EdgeAlpha', 0);
				elseif (iii == 2)
%					aaa = bar(log10(abs(bardata)));
					aaa = bar(abs(bardata));
					XX = get(aaa, 'XData');
					set(aaa, 'FaceColor', 'r', 'EdgeColor', 'r', 'XData', XX, 'FaceAlpha', 1,'LineWidth', 2, 'BarWidth', 0.25, 'EdgeAlpha', 0);
				elseif (iii == 3)
%					aaa = bar(log10(abs(bardata)));
					aaa = bar(abs(bardata));
					XX = get(aaa, 'XData');
					set(aaa, 'FaceColor', 'g', 'EdgeColor', 'g', 'XData', XX + 1/4, 'FaceAlpha', 1, 'LineWidth', 2, 'BarWidth', 0.25,'EdgeAlpha', 0);
				end%if iii
				hold on;
				xlim([1 - 0.5,Nspecial+ 0.5]);
%				ylim([0,5]);
%				ylabel('$\mathrm{log}_{10}|g_i \cdot \partial Y_i \slash \partial Z| [\mathrm{kJ/kg}]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
				ylabel('$|g_i \cdot \partial Y_i \slash \partial Z| [\mathrm{MJ/kg}]$', 'Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
				set(gca, 'FontSize', fs, 'FontName', 'Times', 'Xticklabel', barlabels, 'Xtick', 1:Nspecial);
			end%(bar_graph)
		end%iii = 1:3

		if (kase == 1)
			print -depsc BarGraphLean.eps
		elseif (kase == 2)
			print -depsc BarGraphStoich.eps
		elseif (kase == 3)
			print -depsc BarGraphRich.eps
		end

	end%kase
end%plotSpeciesSensitivity()
