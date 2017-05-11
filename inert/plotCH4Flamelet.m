function[] = plotCH4Flamelet(data)
	data = data';
	Z = data(:,1);
	T = data(:,2);
	Y_CH4 = data(:,18);
	Y_H2O = data(:,31);
	Y_O2 = data(:,5);
	Z_st = 0.0551538;
	x = Z./(Z + Z_st);
	hh = figure();
	set(hh, 'Position', [0 0 650 450]);
	[ax, h1, h2] = plotyy(x, T, x, [Y_CH4, Y_O2, Y_H2O]);
	hold on;
	plot([0.5 0.5], [0 2500],'k--','LineWidth',2)
%	ylim(ax(1), [300, 2400]);
	set(h1, 'LineWidth', 2, 'LineStyle', '-','Color','k');
	set(h2(1), 'LineWidth', 2, 'LineStyle', '--','color', 'b');
	set(h2(2), 'LineWidth', 2, 'LineStyle', '-.','color','r');
	set(h2(3), 'LineWidth', 2, 'LineStyle', ':','color','m');
	xlabel('$Z \slash (Z + Z_{st}) [-]$', 'FontSize', 14, 'FontName', 'Times', 'Interpreter','LaTeX');
	ylabel(ax(1), '$T [K]$','FontSize', 14, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	ylabel(ax(2), '$Y_i [-]$','FontSize', 14, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	set(ax(1), 'FontSize', 14, 'FontName', 'Times','YColor', 'k');
	set(ax(2), 'FontSize', 14, 'FontName', 'Times', 'YColor', 'k');
end
