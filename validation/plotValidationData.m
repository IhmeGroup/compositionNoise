function[] = plotValidationData()
	load('validationData.mat');

	close all;
	h5a = figure();
	h5b = figure();
	h6a = figure();
	h6b = figure();

	figure(h5a);
	plot(OMEGA, abs(TRANS(1,:,1)), 'b-', 'LineWidth', 2); 
	hold on;
	plot(OMEGA, abs(TRANS(2,:,1)), 'r--', 'LineWidth', 2);
	plot(OMEGA, abs(TRANS(3,:,1)), 'k:', 'LineWidth', 2);

	figure(h5b);
	phase1 = PHASE(1,:,1);
	phase2 = PHASE(2,:,1);
	phase3 = PHASE(3,:,1);
%	This unwraps the phase functions so there aren't discontinuous jumps every 2*pi
	phase1 = unwrap(phase1);
	phase2 = unwrap(phase2);
	phase3 = unwrap(phase3);
	plot(OMEGA, phase1, 'b-', 'LineWidth', 2); 
	hold on;
	plot(OMEGA, phase2, 'r--', 'LineWidth', 2);
	plot(OMEGA, phase3, 'k:', 'LineWidth', 2);

	figure(h6a);
	plot(OMEGA, abs(TRANS(1,:,2)), 'b-', 'LineWidth', 2); 
	hold on;
	plot(OMEGA, abs(TRANS(2,:,2)), 'r--', 'LineWidth', 2);
	plot(OMEGA, abs(TRANS(3,:,2)), 'k:', 'LineWidth', 2); 

	figure(h6b);
	phase1 = PHASE(1,:,2);
	phase2 = PHASE(2,:,2);
	phase3 = PHASE(3,:,2);
%	This unwraps the phase functions so there aren't discontinuous jumps every 2*pi
	phase1 = unwrap(phase1);
	phase2 = unwrap(phase2);
	phase3 = unwrap(phase3);

	phase1 = 0*phase1;
	plot(OMEGA, phase1, 'b-', 'LineWidth', 2); 
	hold on;
	plot(OMEGA, phase2, 'r--', 'LineWidth', 2); 
	plot(OMEGA, phase3, 'k:', 'LineWidth', 2);

	figure(h5a);
	hold on;
	red = load('../data/DuranMoreauData/fig5ared.mat');
	plot(red.fig5ared(:,1), red.fig5ared(:,2), 'ro');
	black = load('../data/DuranMoreauData/fig5ablack.mat');
	plot(black.fig5ablack(:,1), black.fig5ablack(:,2), 'ko');
	xlabel('$\Omega$', 'Interpreter', 'laTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('Modulus', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName', 'Times');
	xlim([0 2]);
	ylim([0, 0.6]);
	grid on;

	figure(h5b);
	hold on;
	red = load('../data/DuranMoreauData/fig5bred.mat');
	plot(red.fig5bred(:,1), red.fig5bred(:,2), 'ro');
	black = load('../data/DuranMoreauData/fig5bblack.mat');
	plot(black.fig5bblack(:,1), black.fig5bblack(:,2), 'ko');
	xlabel('$\Omega$', 'Interpreter', 'laTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('Phase', 'FontSize', 14, 'FontName', 'Times');
	xlim([0, 2]);
	ylim([-11/2*pi pi/2]);
	set(gca, 'Xtick', 0:0.4:2.0, 'YTick', -6*pi:pi/2:pi);
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	set(gca, 'YtickLabel', {'$-11\pi \slash 2$', '$-5\pi$','$-9\pi \slash 2$', '$-4\pi$', '$-7\pi \slash 2$', '$-3\pi$','$-5\pi \slash 2$', '$-2\pi$', '$-3\pi \slash 2$', '$-\pi$', '$-1 \pi \slash 2$', '$0$', '$\pi \slash 2$'});
	set(gca, 'FontSize', 14, 'FontName', 'Times');
	grid on;

	figure(h6a);
	hold on;
	red = load('../data/DuranMoreauData/fig6ared.mat');
	plot(red.fig6ared(:,1), red.fig6ared(:,2), 'ro');
	black = load('../data/DuranMoreauData/fig6ablack.mat');
	plot(black.fig6ablack(:,1), black.fig6ablack(:,2), 'ko');
	xlabel('$\Omega$', 'Interpreter', 'laTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('Modulus', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName', 'Times');
	xlim([0, 2]);
	ylim([0,0.9]);
	set(gca, 'XTick', 0:0.4:2.0, 'YTick', 0:0.1:0.9);
	grid on;

	figure(h6b);
	hold on;
	red = load('../data/DuranMoreauData/fig6bred.mat');
	plot(red.fig6bred(:,1), red.fig6bred(:,2), 'ro');
	black = load('../data/DuranMoreauData/fig6bblack.mat');
	plot(black.fig6bblack(:,1), black.fig6bblack(:,2), 'ko');
	xlabel('$\Omega$', 'Interpreter', 'laTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('Phase', 'FontSize', 14, 'FontName', 'Times');
	xlim([0, 2]);
	ylim([-11*pi/2, 3*pi/2])
	set(gca, 'Xtick', 0:0.4:2.0, 'YTick', -11*pi/2:pi/2:3*pi/2);
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	set(gca, 'Ytick', -7*pi:pi:pi, 'YtickLabel', {'$-11\pi \slash 2$', '$-5\pi$', '$-9\pi \slash 2$', '$-4\pi$', '$-7\pi \slash 2$','$-3\pi$', '$-5\pi \slash 2$', '$-2\pi$', '$-3\pi \slash 2$', '$-1\pi$', '$$-1\pi \slash 2$', '$0$', '$1\pi \slash 2$',  '$\pi$', '$3 \pi \slash 2$'});
	set(gca, 'FontSize', 14, 'FontName', 'Times');
	grid on;
end%testOne

