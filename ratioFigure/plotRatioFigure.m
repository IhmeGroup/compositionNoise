function[] = plotRatioFigure();
	load('../EntropicResponseFigure/entropicResponseData.mat');
	TRANS_E = TRANS;
	PHASE_E = PHASE;
	load('../AcousticResponseFigure/acousticResponseData.mat');
	TRANS_A = TRANS;
	PHASE_A = PHASE;
	load('../CompositionResponseFigure/compositionResponseData.mat');
	TRANS_C = TRANS;
	PHASE_C = PHASE;

	close all;
	h5a = figure();
	set(h5a, 'Position', 1.5*[0 0 600, 500]);

	h = 0.225;%figure height
	w = 0.35;%figure width
	lw = 4;%line width for plotting

%	Subplot 1
	subplot('Position', [0.10 0.725 w h]);
	plot(OMEGA, abs(TRANS_C(1,:,1))./abs(TRANS_A(1,:,1)), 'b-', 'LineWidth', lw); 
	hold on;
	plot(OMEGA, abs(TRANS_C(2,:,1))./abs(TRANS_A(1,:,1)), 'r--', 'LineWidth', lw);
	plot(OMEGA, abs(TRANS_C(3,:,1))./abs(TRANS_A(1,:,1)), 'k:', 'LineWidth', lw);
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\frac{|\pi_b^+ \slash \xi_a|}{|\pi_b^+ \slash \pi_a^+|}$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	title('$(a)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	set(gca, 'XTick', 0:0.4:2.0);
	xlim([0 2]);
	grid on;
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(1,1,1))./abs(TRANS_A(1,1,1)), abs(TRANS_C(1,1,1))./abs(TRANS_A(1,1,1))], 'b--', 'LineWidth', lw); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(2,1,1))./abs(TRANS_A(2,1,1)), abs(TRANS_C(2,1,1))./abs(TRANS_A(2,1,1))], 'r--', 'LineWidth', lw); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(3,1,1))./abs(TRANS_A(3,1,1)), abs(TRANS_C(3,1,1))./abs(TRANS_A(3,1,1))], 'k--', 'LineWidth', lw); 

%	Subplot 2
	subplot('Position', [0.55 0.725 w h]);
	plot(OMEGA, abs(TRANS_C(1,:,1))./abs(TRANS_E(1,:,1)), 'b-', 'LineWidth', lw); 
	hold on;
	plot(OMEGA, abs(TRANS_C(2,:,1))./abs(TRANS_E(1,:,1)), 'r--', 'LineWidth', lw);
	plot(OMEGA, abs(TRANS_C(3,:,1))./abs(TRANS_E(1,:,1)), 'k:', 'LineWidth', lw);
	hold on;
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\frac{|\pi_b^+ \slash \xi_a|}{|\pi_b^+ \slash \sigma_a^+|}$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	xlim([0, 2]);
	set(gca, 'XTick', 0:0.4:2.0);
	grid on;
	title('$(b)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(1,1,1))./abs(TRANS_E(1,1,1)), abs(TRANS_C(1,1,1))./abs(TRANS_E(1,1,1))], 'b--', 'LineWidth', lw); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(2,1,1))./abs(TRANS_E(2,1,1)), abs(TRANS_C(2,1,1))./abs(TRANS_E(2,1,1))], 'r--', 'LineWidth', lw); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(3,1,1))./abs(TRANS_E(3,1,1)), abs(TRANS_C(3,1,1))./abs(TRANS_E(3,1,1))], 'k--', 'LineWidth', lw); 

	print -depsc Ratios.eps

end%function
