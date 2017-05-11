function[] = plotValidationData()
	load('compositionResponseData.mat');

	close all;
	h5a = figure();
	set(h5a, 'Position', 1.5*[0 0 600, 500]);

	h = 0.225;%figure height
	w = 0.35;%figure width
	lw = 4;%line width for plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,1)
	plot(OMEGA, abs(TRANS(1,:,1)), 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, abs(TRANS(2,:,1)), 'r', 'LineWidth', 3);
	plot(OMEGA, abs(TRANS(3,:,1)), 'k', 'LineWidth', 4);
	%xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$|\pi_b^+ \slash \xi_a|$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	title('$(a)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	set(gca, 'XTick', 0:0.4:2.0);
	xlim([0 2]);
	grid on;
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(1,1,1)), abs(TRANS(1,1,1))], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(2,1,1)), abs(TRANS(2,1,1))], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(3,1,1)), abs(TRANS(3,1,1))], 'k:', 'LineWidth', 4); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,2)
	plot(OMEGA, abs(TRANS(1,:,2)), 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, abs(TRANS(2,:,2)), 'r', 'LineWidth', 3);
	plot(OMEGA, abs(TRANS(3,:,2)), 'k', 'LineWidth', 4); 
	hold on;
	%xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$|\pi_b^- \slash \xi_a|$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	xlim([0, 2]);
	set(gca, 'XTick', 0:0.4:2.0);
	grid on;
	title('$(b)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(1,1,2)), abs(TRANS(1,1,2))], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(2,1,2)), abs(TRANS(2,1,2))], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(3,1,2)), abs(TRANS(3,1,2))], 'k:', 'LineWidth', 4); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,3)
	plot(OMEGA, abs(TRANS(1,:,5)), 'b', 'LineWidth', 1); 
	hold on;
	plot(OMEGA, abs(TRANS(2,:,5)), 'r', 'LineWidth', 3);
	plot(OMEGA, abs(TRANS(3,:,5)), 'k', 'LineWidth', 4); 
	hold on;
	%xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$|\pi_a^- \slash \xi_a|$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter','LaTeX');
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	xlim([0, 2]);
	set(gca, 'XTick', 0:0.4:2.0);
	grid on;
	title('$(c)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(1,1,5)), abs(TRANS(1,1,5))], 'b:', 'LineWidth', 1); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(2,1,5)), abs(TRANS(2,1,5))], 'r:', 'LineWidth', 3); 
	plot([OMEGA(1), OMEGA(length(OMEGA))], [abs(TRANS(3,1,5)), abs(TRANS(3,1,5))], 'k:', 'LineWidth', 4); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,4)
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
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\angle \pi_b^+ \slash \xi_a$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	xlim([0, 2]);
	set(gca, 'Xtick', 0:0.4:2.0);
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	plot(OMEGA,pi*ones(size(OMEGA)), 'b:', 'LineWidth', 1); 
	plot(OMEGA,pi*ones(size(OMEGA)), 'r:', 'LineWidth', 3); 
	plot(OMEGA,pi*ones(size(OMEGA)), 'k:', 'LineWidth', 4); 
	ylim([-7*pi pi]);
	set(gca, 'YTick', -7*pi:pi:pi);
	set(gca, 'YtickLabel', {'$-7\pi$', '$-6\pi$', '$-5\pi$', '$-4\pi$', '$-3\pi$', '$-2\pi$', '$-\pi$', '0', '$\pi$'});
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	grid on;
	title('$(d)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,5)
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
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\angle \pi_b^- \slash \xi_a$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	xlim([0, 2]);
	set(gca, 'Xtick', 0:0.4:2.0);
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	ylim([-7*pi, 0])
	set(gca, 'YTick', -7*pi:pi:0);
	set(gca, 'YtickLabel', {'$-7\pi$', '$-6\pi$', '$-5\pi$', '$-4\pi$', '$-3\pi$', '$-2\pi$', '$-\pi$', '0'});
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	grid on;
	title('$(e)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Figure 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(3,3,6)
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
	xlabel('$He$', 'Interpreter', 'laTeX', 'FontSize', 18, 'FontName', 'Times');
	ylabel('$\angle \pi_a^- \slash \xi_a$', 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	xlim([0, 2]);
	set(gca, 'Xtick', 0:0.4:2.0);
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	ylim([-pi, 0])
	set(gca, 'YTick', -pi:(pi/2):0);
	set(gca, 'YtickLabel', {'$-\pi$','$-\pi/2$', '0'});
	set(gca, 'FontSize', 18, 'FontName', 'Times');
	grid on;
	title('$(f)$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'Normal','Interpreter', 'LaTeX');

	print -depsc CompositionResponse.eps
end
