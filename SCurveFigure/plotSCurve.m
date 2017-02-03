function[] = plotSCurve(fuel)
%	Plots the Scurve.
	close all;
	hh = figure();
	set(hh, 'Position', [0 0 568 410]*1.1);

%	Load and plot Scurve data
	if (fuel == 1)%C12H26
		data = load('Scurve_C12H26.txt');
	elseif (fuel == 2)%CH4
		data = load('Scurve_CH4.txt');
	elseif (fuel == 3)%H2
		data = load('Scurve_H2.txt');
	end
	data = sortrows(data,4);
	semilogx(data(:,1), data(:,4), 'k-', 'LineWidth',2);
	hold on;

	xlabel('$\chi_{st}$ [s$^{-1}$]', 'Interpreter', 'LaTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$T_{max}$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 18, 'FontName', 'Times');
	xlim([1E-1, 1E2]);
	ylim([600,2300]);
	set(gca, 'FontSize', 18, 'FontName', 'Times');

%	Reload H2 data and plot comparison point
%	data = load('Scurve_C12H26.txt');
	if (fuel == 1)%C12H26
		data = load('Scurve_C12H26.txt');
		plot(data(36,1), data(36,4), 'kx', 'MarkerSize', 18,'LineWidth', 4);
	elseif (fuel == 2)%CH4
%		data = load('Scurve_CH4.txt');
		plot(data(115,1), data(115,4), 'kx', 'MarkerSize', 18,'LineWidth', 4);
		
	elseif (fuel == 3)%H2
		plot(data(73,1), data(73,4), 'kx', 'MarkerSize', 18,'LineWidth', 4);
	end
	print -depsc SCurve.eps
end%plotSCurves()
