function[] = plotSCurve()
%	Plots the Scurve.
	close all;
	hh = figure();
	set(hh, 'Position', [0 0 568 410]*1.1);

%	Load and plot C12H26 data
	data = load('Scurve_C12H26.txt');
	data = sortrows(data,4);
	semilogx(data(:,1), data(:,4), 'k-', 'LineWidth',2);
	hold on;

	xlabel('$\chi_{st}$ [s$^{-1}$]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$T_{max}$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlim([1E-1, 1E2]);
	ylim([600,2300]);
	set(gca, 'FontSize', 14, 'FontName', 'Times');

%	Reload H2 data and plot comparison point
	data = load('Scurve_C12H26.txt');
	plot(data(36,1), data(36,4), 'kx', 'MarkerSize', 14,'LineWidth', 4);
	print -depsc SCurve.eps
end%plotSCurves()
