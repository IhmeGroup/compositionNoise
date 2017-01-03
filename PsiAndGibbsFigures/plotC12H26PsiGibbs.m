function[] = plotC12H26PsiGibbs()
%	Generates 2D plots of Psi and g
	close all;
	addpath('../core');
	addpath('../data');
	global fuel data;
	fuel = 1;%dodecane
	data = loadFuelData(fuel);
	data = data';
	Z = data(:,1);
	T = data(:,2);
	Y_C12H26 = data(:,24);
	Y_H2O = data(:,8);
	Y_O2 = data(:,10);
	Z_st = 0.0627964;
	x = Z./(Z + Z_st);

	p0 = 1E6;
	psi = zeros(length(Z),1);
	g = zeros(length(Z), 1);
	data = data';
	for i = 1:length(Z)
		[psi(i), g(i)] = returnPsi(T(i), p0, Z(i));
	end

	hh = figure();
	set(hh, 'Position', [0 0 560 420]);
	[ax, h1, h2] = plotyy(x, psi, x, g/1E10);
	set(h1, 'LineWidth', 2, 'LineStyle', '-','Color','k');
	set(h2(1), 'LineWidth', 2, 'LineStyle', '--','color', 'b');
	xlabel('$Z \slash (Z + Z_{st}) [-]$', 'FontSize', 14, 'FontName', 'Times', 'Interpreter','LaTeX');
	ylabel(ax(1), '$\Psi [-]$','FontSize', 14, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	ylabel(ax(2), '$g [\textrm{MJ/kg}]$','FontSize', 14, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	set(ax(1), 'FontSize', 14, 'FontName', 'Times','YColor', 'k');
	set(ax(2), 'FontSize', 14, 'FontName', 'Times', 'YColor', 'k');
	set(gca, 'position', [0.16 0.15 0.7 0.75]);
	print -depsc C12H26PsiGibbs.eps
end
