function[] = plotSpeciesSensitivity()
	close all;
%	Flags
	mixture_stuff = true;
	bar_graph = false;

%	Parameters
	NMach = 201;
	T_a = 300;
	p_a = 1.0E6;

%	Load flamelet data
%	data = load('./lowStrain/lowStrain.CH4');
	data = load('./lowStrain/lowStrain.C12H26');
	Z_st = 0.0627964;
	[~, Npts] = size(data);
	[Nspecies, species, a, A, MW] = speciesPropsC12H26();

%	allocate
	PSI = zeros(NMach,1);
	MB = zeros(NMach,1);
	G = zeros(NMach,1);
	GI = zeros(NMach, Nspecies);
	CP = zeros(NMach, 1);
	CPI = zeros(NMach, Nspecies);
	PSII  = zeros(NMach, Nspecies);

%	indices = [14, 26, 81];
%	indices= [9, 26, 99];
	indices = [17, 46, 74];
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

		if (mixture_stuff)
			if (kase == 1)
				hh = figure();
				set(hh, 'Position', [0 0 650 450]);
				plot(MB, PSI, 'b-', 'LineWidth', 2);
				hold on;
				plot(MB, G, 'b--', 'LineWidth', 2);
%				[ax, h1, h2] = plotyy(MB, [PSI], MB, G);
%				set(h1(1), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
%				set(h1(2), 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
%				set(h2, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
%				get(ax(1));
%				hold on;
%				get(ax(2));
%				hold on;
%				ylim(ax(1), [-50, 0]);
%				ylim(ax(2), [1.2, 1.7]);
%				set(ax(2),'FontSize', 14, 'FontName', 'Times', 'YColor', 'k', 'YTick', 1.2:0.1:1.7);
			elseif (kase == 2)
				figure(hh);
				plot(MB, PSI, 'r-','LineWidth',2);
				hold on;
				plot(MB, G, 'r--','LineWidth',2);
%				[ax, h1, h2] = plotyy(ax, MB, [PSI], MB, G);
%				[ax, h1] = plot(MB, [PSI, G]);
%				set(h1(1), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
%				set(h1(2), 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--');
%				set(h2, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--');
%				ylim(ax(1), [-50, 0]);
%				ylim(ax(2), [1.2, 1.7]);
%				set(ax(2),'FontSize', 14, 'FontName', 'Times', 'YColor', 'k', 'YTick', [-1E6,1E6]);
			elseif (kase == 3)
				figure(hh);
				plot(MB, PSI, 'm-','LineWidth',2);
				hold on;
				plot(MB, G, 'm--','LineWidth', 2);
				set(gca,'FontSize', 14, 'FontName', 'Times');
%				[ax, h1, h2] = plotyy(ax, MB, [PSI], MB, G);
%				[ax, h1] = plot(MB, [PSI, G]);
%				set(h1(1), 'LineWidth', 2, 'Color', 'r', 'LineStyle', ':');
%				set(h1(2), 'LineWidth', 2, 'Color', 'b', 'LineStyle', ':');
%				set(h2, 'LineWidth', 2, 'Color', 'b', 'LineStyle', ':');
%				ylim(ax(1), [-50, 0]);
%				ylim(ax(2), [1.2, 1.7]);
%				set(ax(2),'FontSize', 14, 'FontName', 'Times', 'YColor', 'k', 'YTick', [-1E6,1E6]);
			end%kase
			%et(ax(1),'FontSize', 14, 'FontName', 'Times', 'YColor', 'k', 'Ytick', -50:10:0);
			if (kase == 1) 
				xlabel('$M_b [-]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
%				ylabel(ax(1), '$\psi_b [-]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
%				ylabel(ax(2), '$g [\mathrm{MJ/kg}]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
				ylabel('$\psi_b [-], g [\mathrm{MJ/kg}]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
%				ylabel(ax(2), '$c_p [ \mathrm{kJ/kg} ]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
			end%kase == 3
		end%(mixture_stuff)


		if (bar_graph)
			ll = figure();
			set(ll, 'Position', [0 0 650 450]);
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

%				I = [16, 03, 30, 29, 17, 05, 04, 10];
				I = [22, 08, 24, 06, 13, 02, 01];
				Nspecial = length(I');
				bardata = zeros(Nspecial,1);
				barlabels = {};
				for i = 1:Nspecial 
					bardata(i) = psii(I(i));
					barlabels(i) = species(I(i));
				end
				if (iii == 1)
%					aaa = bar(log10(abs(bardata)));
					aaa = bar(abs(bardata));
					XX = get(aaa, 'XData');
					set(aaa, 'FaceColor', 'b', 'EdgeColor', 'b', 'XData', XX + -1/4, 'FaceAlpha', 1, 'LineWidth', 2, 'BarWidth', 0.25, 'EdgeAlpha', 0);
				elseif (iii == 2)
%					aaa = bar(log10(abs(bardata)));
					aaa = bar(abs(bardata));
					XX = get(aaa, 'XData');
					set(aaa, 'FaceColor', 'r', 'EdgeColor', 'r', 'XData', XX, 'FaceAlpha', 1,'LineWidth', 2, 'BarWidth', 0.25, 'EdgeAlpha', 0);
				elseif (iii == 3)
%					aaa = bar(log10(abs(bardata)));
					aaa = bar(abs(bardata));
					XX = get(aaa, 'XData');
					set(aaa, 'FaceColor', 'm', 'EdgeColor', 'm', 'XData', XX + 1/4, 'FaceAlpha', 1, 'LineWidth', 2, 'BarWidth', 0.25,'EdgeAlpha', 0);
				end%if iii
				hold on;
				xlim([1 - 0.5,Nspecial+ 0.5]);
%				ylim([0,5]);
%				ylabel('$\mathrm{log}_{10}|g_i \cdot \partial Y_i \slash \partial Z| [\mathrm{kJ/kg}]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
				ylabel('$|g_i \cdot \partial Y_i \slash \partial Z| [\mathrm{kJ/kg}]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
				set(gca, 'FontSize', 14, 'FontName', 'Times', 'Xticklabel', barlabels, 'Xtick', 1:Nspecial);
			end%(bar_graph)
		end%iii = 1:3
	end%kase
end%plotSpeciesSensitivity()
