function[] = plotRatioFigure();
	load('../EntropicResponseFigure/entropicResponseData.mat');
	TRANS_E = TRANS;
	PHASE_E = PHASE;
    clear TRANS  PHASE; 
	load('../AcousticResponseFigure/acousticResponseData.mat');
	TRANS_A = TRANS;
	PHASE_A = PHASE;
        clear TRANS  PHASE; 
	load('../CompositionResponseFigure/compositionResponseData.mat');
	TRANS_C = TRANS;
	PHASE_C = PHASE;
        clear TRANS  PHASE; 

	close all;
	h5a = figure();
	set(h5a, 'Position', 1.5*[0 0 600, 500]);

	h = 0.225;%figure height
	w = 0.35;%figure width
	lw = 4;%line width for plotting

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot('Position', [0.10 0.725 w h]);
	plot(OMEGA, abs(TRANS_C(1,:,1) )./abs(TRANS_A(1,:,1) ), 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, abs(TRANS_C(2,:,1) + TRANS_C(2,:,2))./abs(TRANS_A(2,:,1) + TRANS_A(2,:,2)), 'r', 'LineWidth', 3);
	plot(OMEGA, abs(TRANS_C(3,:,1) + TRANS_C(3,:,2))./abs(TRANS_A(3,:,1) + TRANS_A(3,:,2)), 'k', 'LineWidth', 4);
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\frac{|(\pi_b^++\pi_b^-) \slash \xi_a|}{|(\pi_b^++\pi_b^-) \slash \pi_a^+|}$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	title('$(a)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	set(gca, 'XTick', 0:0.4:2.0);
	xlim([0 2]);
	grid on;
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(1,1,1) )./abs(TRANS_A(1,1,1) ), abs(TRANS_C(1,1,1) )./abs(TRANS_A(1,1,1) )], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(2,1,1) + TRANS_C(2,1,2))./abs(TRANS_A(2,1,1) + TRANS_A(2,1,2)), abs(TRANS_C(2,1,1) + TRANS_C(2,1,2))./abs(TRANS_A(2,1,1) + TRANS_A(2,1,2))], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(3,1,1) + TRANS_C(3,1,2))./abs(TRANS_A(3,1,1) + TRANS_A(3,1,2)), abs(TRANS_C(3,1,1) + TRANS_C(3,1,2))./abs(TRANS_A(3,1,1) + TRANS_A(3,1,2))], 'k:', 'LineWidth', 4); 

%	Subplot 2
	subplot('Position', [0.55 0.725 w h]);
	plot(OMEGA, abs(TRANS_C(1,:,1) + TRANS_C(1,:,2))./abs(TRANS_E(1,:,1) + TRANS_E(1,:,2)), 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, abs(TRANS_C(2,:,1) + TRANS_C(2,:,2))./abs(TRANS_E(2,:,1) + TRANS_E(2,:,2)), 'r', 'LineWidth', 3);
	plot(OMEGA, abs(TRANS_C(3,:,1) + TRANS_C(3,:,2))./abs(TRANS_E(3,:,1) + TRANS_E(3,:,2)), 'k', 'LineWidth', 4);
	hold on;
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\frac{|(\pi_b^++\pi_b^-) \slash \xi_a|}{|(\pi_b^++\pi_b^-) \slash \sigma_a^+|}$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	xlim([0, 2]);
	set(gca, 'XTick', 0:0.4:2.0);
	grid on;
	title('$(b)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(1,1,1) )./abs(TRANS_E(1,1,1) ), abs(TRANS_C(1,1,1) )./abs(TRANS_E(1,1,1) )], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(2,1,1) + TRANS_C(2,1,2))./abs(TRANS_E(2,1,1) + TRANS_E(2,1,2)), abs(TRANS_C(2,1,1) + TRANS_C(2,1,2))./abs(TRANS_E(2,1,1) + TRANS_E(2,1,2))], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS_C(3,1,1) + TRANS_C(3,1,2))./abs(TRANS_E(3,1,1) + TRANS_E( 3,1,2)), abs(TRANS_C(3,1,1) + TRANS_C(3,1,2))./abs(TRANS_E(3,1,1) + TRANS_E(3,1,2))], 'k:', 'LineWidth', 4); 

	print -depsc Ratios.eps

end%function
