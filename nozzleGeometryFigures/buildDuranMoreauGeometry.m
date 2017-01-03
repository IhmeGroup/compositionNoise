function[] = buildDuranMoreauGeometry()
%	This is a function that builds the spatial coordinate vs. cross-sectional area and spatial coordinate vs Mach number figures for the geometry study at the values beta = -1, 0, and 1
	close all;
	h4a = figure();
	h4b = figure();

%	Pre-compute for convenience
	gamma = 1.4;
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;

	for test = 1:3
		if (test == 1)	M_a = 0.29; M_b = 0.88;
		elseif (test == 2) M_a = 0.29; M_b = 1.02;
		elseif (test == 3) M_a = 0.29; M_b = 1.5;
		end

%		Compute the Mach number and cross sectional area
		N = 1000;
		etahat = linspace(gp1o2*M_a.^2./(1 + gm1o2*M_a.^2), gp1o2*M_b.^2./(1+gm1o2*M_b.^2), N);
		M = sqrt(2/gp1.*etahat./(1 - gm1/gp1.*etahat));
		A = gp1o2.^(-gp1o2/gm1).*(1 + gm1o2.*M.*M).^(gp1o2/gm1)./M;
		eta = sqrt(etahat);
		eta = (eta - min(eta));
		eta  = eta./max(eta);

%		plot the cross sectional area vs. eta
		figure(h4a);
		if (test == 1)		plot(eta, A, 'b-', 'LineWidth', 2);	hold on;
		elseif (test == 2)	plot(eta, A, 'r--', 'LineWidth', 2);
		elseif (test == 3)	plot(eta, A, 'k:', 'LineWidth', 2);
		end

		figure(h4b);
		if (test == 1)		plot(eta, M, 'b-', 'LineWidth', 2);	hold on;
		elseif (test == 2)	plot(eta, M, 'r--', 'LineWidth', 2);
		elseif (test == 3)	plot(eta, M, 'k:', 'LineWidth', 2);
		end
	end

%	Format the eta vs. A figure
	figure(h4a);
	ylim([-0,2.5]);
	xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$A(\eta) \slash A^*$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'Xtick', [0:0.2:1.0], 'YTick', [0:0.2:2.4]);
	grid on;

%	Plot the profile
	figure(h4b);
	ylim([0, 1.8]);
	xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$\bar{M}(\eta)$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'Xtick', [0:0.2:1.0], 'YTick', [0:0.2:1.8]);
	grid on;

end%buildNozzleGeometry
