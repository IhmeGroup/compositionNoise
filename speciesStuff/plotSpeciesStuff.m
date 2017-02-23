function[] = plotSpeciesStuff()
	[Nspecies, species, a, A, MW] = speciesPropsCH4();
	close all;

	lw = 4;
	fs = 22;


	NMach = 371;
	Tvec = [1320, 2030, 1710];
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
	ymax = 0;

%		FUEL OXYGEN NITROGEN WATER CO2 RADICALS		
	I = [16, 03, 30, 29, 17, 05];
%	gg = figure();
%	set(gg, 'Position', [0 0 650 450]);
	figure();
	h1 = plot(TB, [GI(:,I(1)), GI(:,I(2)), GI(:,I(3)), GI(:,I(4)), GI(:,I(5)), GI(:,I(6))]);
	hold on;
	hold on;
%	plot([1320, 1320], [ymin,ymax], 'b-', 'LineWidth', 3);
%	plot([2030, 2030], [ymin,ymax], 'r-', 'LIneWidth', 2);
%	plot([1710, 1710], [ymin,ymax], 'm-', 'LineWidth', 1);
%	plot([1320, 1320].*.8333, [ymin,ymax], 'b--', 'LineWidth', 3);
%	plot([2030, 2030].*.8333, [ymin,ymax], 'r--', 'LIneWidth', 2);
%	plot([1710, 1710].*.8333, [ymin,ymax], 'm--', 'LineWidth', 1);
%	plot([1320, 1320].*.5556, [ymin,ymax], 'b:', 'LineWidth', 3);
%	plot([2030, 2030].*.5556, [ymin,ymax], 'r:', 'LIneWidth', 2);
%	plot([1710, 1710].*.5556, [ymin,ymax], 'm:', 'LineWidth', 1);
	set(gca, 'FontSize', 16, 'Fontname', 'Times','YColor', 'k');
	set(h1(1), 'LineWidth', lw, 'Color', 'k', 'LineStyle', '-');
	set(h1(2), 'LineWidth', lw, 'Color', 'r', 'LineStyle', '--');
	set(h1(3), 'LineWidth', lw, 'Color', 'm', 'LineStyle', '-.');
	set(h1(4), 'LineWidth', lw, 'Color', 'b', 'LineStyle', ':');
	set(h1(5), 'LineWidth', lw/4, 'Color', 'c', 'LineStyle', '-');
	set(h1(6), 'LineWidth', lw/4, 'Color', 'r', 'LineStyle', '--');
	xlabel('$T_c [K]$', 'Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
	ylabel('$g_i [MJ/kg]$', 'Interpreter','LaTeX', 'FontSize', fs, 'FontName', 'Times');
	set(gca, 'FontSize', fs, 'FontName', 'Times');
	jj = legend([species(I(1)), species(I(2)), species(I(3)), species(I(4)), species(I(5)), species(I(6))]);
	set(jj,'Location', 'SouthWest', 'FontSize', fs);
	xlim([700, 2100]);
	ylim([ymin, ymax]);

	print -djpeg SpeciesGibbs.jpg
end
