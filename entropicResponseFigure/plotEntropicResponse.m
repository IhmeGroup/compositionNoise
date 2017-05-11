function[] = plotValidationData()
	load('entropicResponseData.mat');

	close all;
	h5a = figure();
	set(h5a, 'Position', 1.5*[0 0 600, 500]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,1);
	plot(OMEGA, abs(TRANS(1,:,1)), 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, abs(TRANS(2,:,1)), 'r', 'LineWidth', 3);
	plot(OMEGA, abs(TRANS(3,:,1)), 'k', 'LineWidth', 4);
	red = load('../data/DuranMoreauData/fig5ared.mat');
    
	plot(red.fig5ared(:,1), red.fig5ared(:,2), 'w-o', 'MarkerEdgeColor', 'r','MarkerSize',10);
	black = load('../data/DuranMoreauData/fig5ablack.mat');
	plot(black.fig5ablack(:,1), black.fig5ablack(:,2), 'ko','MarkerSize',10);
    
	%xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$|\pi_b^+ \slash \sigma_a|$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	title('$(a)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	xlim([0 2]);
	ylim([0, 0.6]);
	grid on;
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(1,1,1)), abs(TRANS(1,1,1))], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(2,1,1)), abs(TRANS(2,1,1))], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(3,1,1)), abs(TRANS(3,1,1))], 'k:', 'LineWidth', 4); 
%	Legend
	ll = legend('$M_b$ = 0.88 In-house', '$M_b$ = 1.02 In-house', '$M_b$ = 1.50 In-house', '$M_b$ = 1.02 Duran \& Moreau', '$M_b$ = 1.50 Duran \& Moreau', '$M_b$ = 0.88 Compact Nozzle', '$M_b$ = 1.02 Compact Nozzle', '$M_b$ = 1.50 Compact Nozzle');
	set(ll, 'Position', [0.35 0.05 0.30 0.3], 'Interpreter', 'LaTeX', 'FontName', 'Times', 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,2);
	plot(OMEGA, abs(TRANS(1,:,2)), 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, abs(TRANS(2,:,2)), 'r', 'LineWidth', 3);
	plot(OMEGA, abs(TRANS(3,:,2)), 'k', 'LineWidth', 4); 
	hold on;
    
	red = load('../data/DuranMoreauData/fig6ared.mat');
	plot(red.fig6ared(:,1), red.fig6ared(:,2), 'w-o', 'MarkerEdgeColor', 'r','MarkerSize',10);
	black = load('../data/DuranMoreauData/fig6ablack.mat');
	plot(black.fig6ablack(:,1), black.fig6ablack(:,2), 'ko','MarkerSize',10);
	%xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$|\pi_b^- \slash \sigma|$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	xlim([0, 2]);
	set(gca, 'Xtick', 0:0.4:2.0);
	ylim([0,0.9]);
	set(gca, 'XTick', 0:0.4:2.0, 'YTick', 0:0.2:1.0);
	grid on;
	title('$(b)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(1,1,2)), abs(TRANS(1,1,2))], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(2,1,2)), abs(TRANS(2,1,2))], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(3,1,2)), abs(TRANS(3,1,2))], 'k:', 'LineWidth', 4); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,3);
	plot(OMEGA, abs(TRANS(1,:,5)), 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, abs(TRANS(2,:,5)), 'r', 'LineWidth', 3);
	plot(OMEGA, abs(TRANS(3,:,5)), 'k', 'LineWidth', 4); 
	hold on;
    
	%red = load('../data/DuranMoreauData/fig6ared.mat');
	%plot(red.fig6ared(:,1), red.fig6ared(:,2), 'w-o', 'MarkerEdgeColor', 'r','MarkerSize',10);
	%black = load('../data/DuranMoreauData/fig6ablack.mat');
	%plot(black.fig6ablack(:,1), black.fig6ablack(:,2), 'ko','MarkerSize',10);
	%xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$|\pi_a^- \slash \sigma|$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	xlim([0, 2]);
	set(gca, 'Xtick', 0:0.4:2.0);
	ylim([0,0.3]);
	set(gca, 'XTick', 0:0.4:2.0, 'YTick', 0:0.1:0.3);
	grid on;
	title('$(c)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(1,1,5)), abs(TRANS(1,1,5))], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(2,1,5)), abs(TRANS(2,1,5))], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(3,1,5)), abs(TRANS(3,1,5))], 'k:', 'LineWidth', 4); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,4);
	phase1 = PHASE(1,:,1);
	phase2 = PHASE(2,:,1);
	phase3 = PHASE(3,:,1);
%	This unwraps the phase functions so there aren't discontinuous jumps every 2*pi
	phase1 = unwrap(phase1);
	phase2 = unwrap(phase2);
	phase3 = unwrap(phase3);
	plot(OMEGA, phase1, 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, phase2, 'r', 'LineWidth', 3);
	plot(OMEGA, phase3, 'k', 'LineWidth', 4);
	red = load('../data/DuranMoreauData/fig5bred.mat');
	plot(red.fig5bred(:,1), red.fig5bred(:,2), 'w-o', 'MarkerEdgeColor', 'r','MarkerSize',10);
	black = load('../data/DuranMoreauData/fig5bblack.mat');
	plot(black.fig5bblack(:,1), black.fig5bblack(:,2), 'ko','MarkerSize',10);
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\angle \pi_b^+ \slash \sigma_a$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	xlim([0, 2]);
	ylim([-11/2*pi pi/2]);
	set(gca, 'Xtick', 0:0.4:2.0);
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'k:', 'LineWidth', 4); 
	set(gca, 'YTick', -6*pi:pi:pi, 'YTicklabel', {'$-6\pi$', '$-5\pi$' '$-4\pi$', '$-3\pi$', '$-2\pi$', '$-\pi$', '0', '$\pi$'});
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	grid on;
	title('$(d)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	subplot(3,3,5);
	phase1 = PHASE(1,:,2);
	phase2 = PHASE(2,:,2);
	phase3 = PHASE(3,:,2);
%	This unwraps the phase functions so there aren't discontinuous jumps every 2*pi
	phase1 = unwrap(phase1);
	phase2 = unwrap(phase2);
	phase3 = unwrap(phase3);
	phase1 = 0*phase1;
	plot(OMEGA, phase1, 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, phase2, 'r', 'LineWidth', 3); 
	plot(OMEGA, phase3, 'k', 'LineWidth', 4);
    
	red = load('../data/DuranMoreauData/fig6bred.mat');
	plot(red.fig6bred(:,1), red.fig6bred(:,2), 'w-o', 'MarkerEdgeColor', 'r','MarkerSize',10);
	black = load('../data/DuranMoreauData/fig6bblack.mat');
	plot(black.fig6bblack(:,1), black.fig6bblack(:,2), 'ko','MarkerSize',10);
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\angle \pi_b^- \slash \sigma_a$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	xlim([0, 2]);
	ylim([-11*pi/2, 3*pi/2])
	set(gca, 'Xtick', 0:0.4:2.0, 'YTick', -11*pi/2:pi/2:3*pi/2);
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	set(gca, 'Ytick', -7*pi:pi:pi, 'Yticklabel', {'$-7\pi$', '$-6\pi$', '$-5\pi$', '$-4\pi$', '$-3\pi$', '$-2\pi$', '$-\pi$', '0', '$\pi$'});
	set(gca, 'FontSize', 18, 'FontName', 'Times');
    
%     plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'b:', 'LineWidth', 1); 
% 	plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'r:', 'LineWidth', 3); 
% 	plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'k:', 'LineWidth', 4); 
    
    plot(OMEGA, 0*ones(size(OMEGA)), 'b:', 'LineWidth', 1); 
	plot(OMEGA, pi*ones(size(OMEGA)), 'r:', 'LineWidth', 3); 
	plot(OMEGA, pi*ones(size(OMEGA)), 'k:', 'LineWidth', 4); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	subplot(3,3,6);
	phase1 = PHASE(1,:,3);
	phase2 = PHASE(2,:,3);
	phase3 = PHASE(3,:,3);
%	This unwraps the phase functions so there aren't discontinuous jumps every 2*pi
	phase1 = unwrap(phase1);
	phase2 = unwrap(phase2);
	phase3 = unwrap(phase3);
	phase1 = 0*phase1;
	plot(OMEGA, phase1, 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, phase2, 'r', 'LineWidth', 3); 
	plot(OMEGA, phase3, 'k', 'LineWidth', 4);
    
	%red = load('../data/DuranMoreauData/fig6bred.mat');
	%plot(red.fig6bred(:,1), red.fig6bred(:,2), 'w-o', 'MarkerEdgeColor', 'r','MarkerSize',10);
	%black = load('../data/DuranMoreauData/fig6bblack.mat');
	%plot(black.fig6bblack(:,1), black.fig6bblack(:,2), 'ko','MarkerSize',10);
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\angle \pi_a^- \slash \sigma_a$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	xlim([0, 2]);
	ylim([0, pi*1.1])
	set(gca, 'Xtick', 0:0.4:2.0, 'YTick', 0:pi/2:pi);
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	set(gca, 'Ytick', 0:pi/2:pi, 'Yticklabel', {'0', '$\pi/2$', '$\pi$'});
	set(gca, 'FontSize', 18, 'FontName', 'Times');
    
%     plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'b:', 'LineWidth', 1); 
% 	plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'r:', 'LineWidth', 3); 
% 	plot([OMEGA(1), OMEGA(length(OMEGA))], [0 0], 'k:', 'LineWidth', 4); 
    
    plot(OMEGA, 0*ones(size(OMEGA)), 'b:', 'LineWidth', 1); 
	plot(OMEGA, pi*ones(size(OMEGA)), 'r:', 'LineWidth', 3); 
	plot(OMEGA, pi*ones(size(OMEGA)), 'k:', 'LineWidth', 4); 
    
	grid on;
	title('$(f)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');

	print -depsc EntropicResponse.eps
end
