function[] = plotSCurves()
	lw = 4;
	fs = 22;
	hh = figure();
	set(hh, 'Position', [0 0 568 410]*1.1);
%	Load and plot H2 data
	data = load('Scurve_H2.txt');
	data = sortrows(data, 4);
	semilogx(data(:,1), data(:,4),'b-','LineWidth',2);
	hold on;

%	Load and plot CH4 data
	data = load('Scurve_CH4.txt');
	data = sortrows(data,4);
	semilogx(data(:,1), data(:,4),'r--','LineWidth',2);

%	Load and plot C12H26 data
	data = load('Scurve_C12H26.txt');
	data = sortrows(data,4);
	semilogx(data(:,1), data(:,4), 'm-.', 'LineWidth',2);

	xlabel('$\chi_{\textrm{st}}$ [s$^{-1}$]', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
	ylabel('$T_{\textrm{st}}$ [K]', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
	xlim([1E-1, 1E3]);
	set(gca, 'FontSize', fs, 'FontName', 'Times');
	aa = legend('$H_2$', '$CH_4$', '$C_{12}H_{26}$');

%	Reload H2 data and plot comparison point
	data = load('Scurve_H2.txt');
	plot(data(4,1), data(4,4), 'bx', 'MarkerSize', fs,'LineWidth', lw);
	data = load('Scurve_CH4.txt');
	plot(data(77,1), data(77,4), 'rx', 'MarkerSize', fs,'LineWidth', lw);
	data = load('Scurve_C12H26.txt');
	plot(data(36,1), data(36,4), 'mx', 'MarkerSize', fs,'LineWidth', lw);

	set(aa, 'FontSize', fs, 'FontName','Times', 'Interpreter','LaTeX');

	plot([1 1], [300,2500], 'k--', 'LineWidth', lw);
	ylim([500,2500]);

	print -depsc SCurves.eps
end%plotSCurves()
