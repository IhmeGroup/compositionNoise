function[eta A M] = validateDuranMoreau()
%	This function is intended to validate the (non-composition) portion of the solver by replicating the data found in
%	Duran and Moreau, Solution of the quasi-one-dimensional linearized Euler equations using flow invariants and the Magnus expansion, JFM vol 723, pp 190-231. 2013
%	The figure number (i.e. 5a, 5b, etc.) correspond to the equivalent figure in Duran and Moreau
%	After generating the equivalent data using O'Brien's solver, results pulled from D&M's paper are then superimposed on the plot at the end of the code
%	This data was taken using grabit in MATLAB and is thus fairly imprecise. Discrepencies between D&M and this solver can be attributed partially to this and partially to numerical error.
%	Varying the ``tightness'' of the throat perturbation parameter epsilon is found to have an appreciable effect on the agreement (small values of epsilon are better), but this also slows down the code
%	Also, a subsonic case is included b/c it seemed like a valuable addition.

	close all;
	h5a = figure();
	h5b = figure();
	h6a = figure();
	h6b = figure();

%	Flags
	use_runway = false;

%	parameters
	gamma = 1.4;
	R_univ = 8.31415;
	Nomega = 31;
	dOmega = 2/(Nomega-1);
	global fuel;
	fuel = 2;

	for test = 1:3
		if (test == 1)		M_a = 0.29; M_b = 0.88; M_c = 0.0;
		elseif (test == 2) 	M_a = 0.29; M_b = 1.02; M_c = 0.0; 
		elseif (test == 3) 	M_a = 0.29; M_b = 1.50; M_c = 0.0; 
		end

		[SPLINES] = buildBaseFlowSplines();
%		Build the base flow
		count = 0;
		for omega = 0:dOmega:2.0
			disp(omega)
			count = count + 1;
			OMEGA(count) = omega;
			if ((exist('subsol', 'var')) && (exist('supsol', 'var')))
				[transfer, subsol, supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, M_c, omega, fuel, true, SPLINES, subsol, supsol);
			else
				[transfer, subsol, supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, M_c, omega, fuel, true, SPLINES);
			end
%			 1         2       3       4  5  6  7  8  9
			TRANS(count,:) = [transfer(1,2), transfer(2,2), transfer(3,2), transfer(4,2), transfer(2,1)];
			PHASE(count) = [atan2(imag(transfer(1,2)), real(transfer(1,2)))];
			PHASE2(count) = [atan2(imag(transfer(2,2)), real(transfer(2,2)))];
			PHASE3(count) = [atan2(imag(transfer(2,1)), real(transfer(2,1)))];
		end%omega
		clear subsol supsol;

		TRANS(:,5)'

		figure(h5a);
		if (test == 1) plot(OMEGA, abs(TRANS(:,1)), 'b-', 'LineWidth', 2); hold on;
		elseif (test == 2) plot(OMEGA, abs(TRANS(:,1)), 'r--', 'LineWidth', 2);
		elseif (test == 3) plot(OMEGA, abs(TRANS(:,1)), 'k:', 'LineWidth', 2); end

		figure(h5b);
		phase = PHASE;
%		This unwraps the phase functions so there aren't discontinuous jumps every 2*pi
		dphase = diff(phase);
		for i = 1:length(dphase)
			if (abs(dphase(i)) > pi)
				phase(i+1:end) = phase(i+1:end) - 2*pi;
				dphase = diff(phase);
			end
		end
		if (test == 1) plot(OMEGA, phase, 'b-', 'LineWidth', 2); hold on;
		elseif (test == 2) plot(OMEGA, phase, 'r--', 'LineWidth', 2);
		elseif (test == 3) plot(OMEGA, phase, 'k:', 'LineWidth', 2); end

		figure(h6a);
		if (test == 1) plot(OMEGA, abs(TRANS(:,2)), 'b-', 'LineWidth', 2); hold on;
		elseif (test == 2) plot(OMEGA, abs(TRANS(:,2)), 'r--', 'LineWidth', 2);
		elseif (test == 3) plot(OMEGA, abs(TRANS(:,2)), 'k:', 'LineWidth', 2); end

		figure(h6b);
		phase = PHASE2;
%		This unwraps the phase functions so there aren't discontinuous jumps every 2*pi
		dphase = diff(phase);
		for i = 1:length(dphase)
			if (abs(dphase(i)) > pi)
				phase(i+1:end) = phase(i+1:end) - 2*pi;
				dphase = diff(phase);
			end
		end
		if (test == 1) phase = 0*phase;end
		if (test == 1) plot(OMEGA, phase, 'b-', 'LineWidth', 2); hold on;
		elseif (test == 2) plot(OMEGA, phase, 'r--', 'LineWidth', 2); 
		elseif (test == 3) plot(OMEGA, phase, 'k:', 'LineWidth', 2); end
	end%for test


	figure(h5a);
	hold on;
	red = load('./DuranMoreauData/fig5ared.mat');
	plot(red.fig5ared(:,1), red.fig5ared(:,2), 'ro');
	black = load('./DuranMoreauData/fig5ablack.mat');
	plot(black.fig5ablack(:,1), black.fig5ablack(:,2), 'ko');
	xlabel('$\Omega$', 'Interpreter', 'laTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('Modulus', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName', 'Times');
	xlim([0 2]);
	ylim([0, 0.6]);
	grid on;

	figure(h5b);
	hold on;
	red = load('./DuranMoreauData/fig5bred.mat');
	plot(red.fig5bred(:,1), red.fig5bred(:,2), 'ro');
	black = load('./DuranMoreauData/fig5bblack.mat');
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
	red = load('./DuranMoreauData/fig6ared.mat');
	plot(red.fig6ared(:,1), red.fig6ared(:,2), 'ro');
	black = load('./DuranMoreauData/fig6ablack.mat');
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
	red = load('./DuranMoreauData/fig6bred.mat');
	plot(red.fig6bred(:,1), red.fig6bred(:,2), 'ro');
	black = load('./DuranMoreauData/fig6bblack.mat');
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

