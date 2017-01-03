function[] = plotHelmholtzSweepData();
	close all;
	load('HelmholtzSweepData.mat');
	printem = true;
	style = {'k-'; 'b-'; 'r-'; 'm-'; 'c-'};
	h1 = figure();
	h2 = figure();
	h3 = figure();
	h4 = figure();
	g1 = figure();
	g2 = figure();
	g3 = figure();
	g4 = figure();

	for i = 1:5
		if i == 1 style = 'k-'; end
		if i == 2 style = 'b-'; end
		if i == 3 style = 'r-'; end
		if i == 4 style = 'm-'; end
		if i == 5 style = 'c-'; end
		eta = cell2mat(ETA(i));
		w_p = cell2mat(WP(i));
		w_m = cell2mat(WM(i));
		w_s = cell2mat(WS(i));
		w_z = cell2mat(WZ(i));


%		This transforms the eta coordinate s.t. eta is always on [0, 1], which makes it easier to compare cases at different Mach numbers
		eta = eta - min(eta);
		eta = eta./max(eta);

		figure(h1);
		plot(eta, abs(w_p), style, 'LineWidth', 2);
		hold on;
		figure(h2);
		plot(eta, abs(w_m), style, 'LineWidth', 2);
		hold on;
		figure(h3);
		plot(eta, abs(w_s), style, 'LineWidth', 2);
		hold on;
		figure(h4);
		plot(eta, abs(w_z), style, 'LineWidth', 2);
		hold on;

		phase_p = atan2(imag(w_p), real(w_p));
		phase_m = atan2(imag(w_m), real(w_m));
		phase_s = atan2(imag(w_s), real(w_s));
		phase_z = atan2(imag(w_z), real(w_z));
		phase_p(1) = 0;
		phase_s = 0*phase_s;

%		This unwraps the function to avoid discontinuous 2\pi jumps in the phase plots.
%		the conditioning on eta for w_m is to avoid the oscillations associated with the singularity

		for j = 2:length(phase_p)
			if (abs(phase_p(j) - phase_p(j-1)) > pi)
				phase_p(j:end) = phase_p(j:end) - 2*pi;
			end
		end

		for j = 2:length(phase_m)
			if (abs(phase_m(j) - phase_m(j-1)) > 0.9*pi)
%				if (eta(j) < 0.99) || (eta(j) > 1.01)
				if (eta(j) < 0.62) || (eta(j) > 0.66)
					phase_m(j:end) = phase_m(j:end) - 2*pi;
				end
			end
		end

		for j = 2:length(phase_s)
			if (abs(phase_s(j) - phase_s(j-1)) > pi)
				phase_s(j:end) = phase_s(j:end) - 2*pi;
			end
		end

		for j = 2:length(phase_z)
			if (abs(phase_z(j) - phase_z(j-1)) > pi)
				phase_z(j:end) = phase_z(j:end) - 2*pi;
			end
		end

		figure(g1);
		plot(eta(2:end), phase_p(2:end), style, 'LineWidth', 2);
		hold on;
		figure(g2);
		plot(eta, phase_m, style, 'LineWidth', 2);
		hold on;
		figure(g3);
		plot(eta, phase_s, style, 'LineWidth', 2);
		hold on;
		figure(g4);
		plot(eta, phase_z, style, 'LineWidth', 2);
		hold on;
	end%for i = 1:%		

	figure(h1);
	xlim([min(eta), max(eta)]);
	set(gca,'FontSize', 14, 'FontName', 'Times');
	xlabel('$\eta$', 'Interpreter', 'LaTeX','FontSize', 14, 'FontName', 'Times');
	ylabel('$|\pi^+(\eta) \slash \xi_a|$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

	figure(h2);
	xlim([min(eta), max(eta)]);
	ylim([0, 1.5]);
	set(gca,'FontSize', 14, 'FontName', 'Times');
	xlabel('$\eta$', 'Interpreter', 'LaTeX','FontSize', 14, 'FontName', 'Times');
	ylabel('$|\pi^-(\eta) \slash xi_a|$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

	figure(h3);
	xlim([min(eta), max(eta)]);
	ylim([0,1.0]);
	set(gca,'FontSize', 14, 'FontName', 'Times');
	xlabel('$\eta$', 'Interpreter', 'LaTeX','FontSize', 14, 'FontName', 'Times');
	ylabel('$|\sigma(\eta) \slash \xi_a|$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

	figure(h4);
	xlim([min(eta), max(eta)]);
	ylim([0,1.1]);
	set(gca,'FontSize', 14, 'FontName', 'Times');
	xlabel('$\eta$', 'Interpreter', 'LaTeX','FontSize', 14, 'FontName', 'Times');
	ylabel('$|\xi(\eta) \slash \xi_a|$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

	figure(g1);
	xlim([min(eta), max(eta)]);
	ylim([-6.5*pi, pi]);
	set(gca, 'YTick', [-6*pi:pi:0], 'YTickLabel', {'-6\pi', '-5\pi', '-4\pi', '-3\pi', '-2\pi', '-1\pi','0\pi', '1\pi'});
	set(gca,'FontSize', 14, 'FontName', 'Times');
	xlabel('$\eta$', 'Interpreter', 'LaTeX','FontSize', 14, 'FontName', 'Times');
	ylabel('$\angle \pi^+(\eta) \slash \xi_a$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

	figure(g2);
	xlim([min(eta), max(eta)]);
	ylim([-7*pi, pi]);
	set(gca, 'YTick', [-7*pi:pi:pi], 'YTickLabel', {'-7\pi', '-6\pi', '-5\pi', '-4\pi', '-3\pi', '-2\pi', '-1\pi','0\pi', '\pi'});
	set(gca,'FontSize', 14, 'FontName', 'Times');
	xlabel('$\eta$', 'Interpreter', 'LaTeX','FontSize', 14, 'FontName', 'Times');
	ylabel('$\angle \pi^-(\eta) \slash \xi_a$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

	figure(g3);
	xlim([min(eta), max(eta)]);
	ylim([-pi, pi]);
	set(gca, 'YTick', [-pi:pi/2:pi], 'YTickLabel', {'$-\pi$', '$-\pi \slash 2$', '0', '$\pi \slash 2$', '$\pi$'});
	set(gca, 'TickLabelInterpreter', 'LaTeX');
	set(gca,'FontSize', 14, 'FontName', 'Times');
	xlabel('$\eta$', 'Interpreter', 'LaTeX','FontSize', 14, 'FontName', 'Times');
	ylabel('$\angle \sigma(\eta) \slash \xi_a$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');

	figure(g4);
	xlim([min(eta), max(eta)]);
	ylim([-6.5*pi, 0.5*pi]);
	set(gca, 'YTick', [-6*pi:pi:0], 'YTickLabel', {'-6\pi', '-5\pi', '-4\pi', '-3\pi', '-2\pi', '-1\pi','0\pi'});
	set(gca,'FontSize', 14, 'FontName', 'Times');
	xlabel('$\eta$', 'Interpreter', 'LaTeX','FontSize', 14, 'FontName', 'Times');
	ylabel('$\angle \xi(\eta) \slash \xi_a$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
end%HelmholtzSweep()
