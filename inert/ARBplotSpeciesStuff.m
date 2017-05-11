function[] = ARBplotSpeciesStuff()
	[Nspecies, species, a, A, MW] = speciesPropsC12H26();
	close all;
	NMach = 371;
	index = [17, 46, 74];
	Tvec = [1293, 2104,  1431];
%	Tvec = [1320, 2030, 1710];
	TB = zeros(NMach,1);
	GI = zeros(NMach, Nspecies);
	T_a = 300.0;
	p_b = 1E6;
	for i = 1:NMach
		T_b = 300 + (i-1)*5;
		TB(i) = T_b;
		
		for j = 1:Nspecies
			Y = zeros(Nspecies,1);
			Y(j) = 1.0;
			[~, ~, ~, GI(i,j)] = returnSpeciesProperties(T_b, p_b, Y, a, A, MW);
		end%j = 1:Nspecies
	end%i = 1:NMach	
%	1 erg/gram is 1E-7 kJ/kg
	GI = GI*1E-10;

	ymin = -40;
	ymax = 10;

%		FUEL OXYGEN NITROGEN WATER CO2 RADICALS		
%	I = [16, 03, 30, 29, 17, 05];
	I = [22, 08, 24, 06, 13, 02, 01];
	gg = figure();
	set(gg, 'Position', [0 0 650 450]);
	h1 = plot(TB, [GI(:,I(1)), GI(:,I(2)), GI(:,I(3)), GI(:,I(4)), GI(:,I(5)), GI(:,I(6))]);
	hold on;
	plot([1320, 1320], [ymin,ymax], 'b-', 'LineWidth', 2);
	plot([2030, 2030], [ymin,ymax], 'r-', 'LIneWidth', 2);
	plot([1710, 1710], [ymin,ymax], 'm-', 'LineWidth', 2);
	plot([1320, 1320].*.8333, [ymin,ymax], 'b--', 'LineWidth', 2);
	plot([2030, 2030].*.8333, [ymin,ymax], 'r--', 'LIneWidth', 2);
	plot([1710, 1710].*.8333, [ymin,ymax], 'm--', 'LineWidth', 2);
	plot([1320, 1320].*.5556, [ymin,ymax], 'b:', 'LineWidth', 2);
	plot([2030, 2030].*.5556, [ymin,ymax], 'r:', 'LIneWidth', 2);
	plot([1710, 1710].*.5556, [ymin,ymax], 'm:', 'LineWidth', 2);
	set(gca, 'FontSize', 14, 'Fontname', 'Times','YColor', 'k');
	set(h1(1), 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
	set(h1(2), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
	set(h1(3), 'LineWidth', 2, 'Color', 'm', 'LineStyle', '-');
	set(h1(4), 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
	set(h1(5), 'LineWidth', 2, 'Color', 'c', 'LineStyle', '-');
	set(h1(6), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '-');
	xlabel('$T_b [K]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$g [MJ/kg]$', 'Interpreter','LaTeX', 'FontSize', 14, 'FontName', 'Times');
	jj = legend([species(I(1)), species(I(2)), species(I(3)), species(I(4)), species(I(5)), species(I(6))]);
	set(jj,'Location', 'SouthWest');
	xlim([700, 2100]);
	ylim([ymin, ymax]);
end
