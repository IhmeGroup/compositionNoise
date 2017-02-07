function[] = plotNozzleData()
	close all;
	load('nozzleData.mat');

	lw = 4;

	figure();
	plot(BETA, abs(TRANS(:, 1, 1)), 'k-', 'LineWidth', lw);
	hold on;
	plot(BETA, abs(TRANS(:, 2, 1)), 'r--', 'LineWidth', lw);
	plot(BETA, abs(TRANS(:, 3, 1)), 'b:', 'LineWidth', lw);
	plot(BETA, abs(TRANS(:, 4, 1)), 'm-.', 'LineWidth', lw);
	plot(BETA, abs(TRANS(:, 5, 1)), 'c-', 'LineWidth', 0.5*lw);
	xlabel('$\beta$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$|\pi_b^+ \slash \xi_a|$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
%	title('$(a)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
%	set(gca, 'XTick', 0:0.4:2.0);
%	xlim([0 2]);
	grid on;
	ll = legend('$He = 0.0$', '$He = 0.5$', '$He = 1.0$', '$He = 1.5$', '$He = 2.0$');
	set(ll, 'Interpreter', 'LaTeX', 'FontSize', 18, 'FontName', 'Times');
	set(ll, 'Location', 'SouthEast');

	print -depsc GeometricSensitivity.eps
end%function
