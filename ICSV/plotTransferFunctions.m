function[] = plotTransferFunctions()
	close all;
	ms = 14;
	for i = 1:3
		if (i == 1)
			load('Acoustic.mat');
			jeonglae = 		[	0.1	1.36535
								0.5	2.00194
								1.0	2.36684
								2.0	2.59787];

		elseif (i == 2)
			load('Entropy.mat');
			jeonglae = 		[	0.1	0.346219	0.360221
								0.5	0.205958	0.243254
								1.0	0.144151	0.155589
								2.0	0.0984836	0.0759971];
		elseif (i == 3)
			load('Composition.mat');
			jeonglae = 		[	0.1	8.1925		8.59796;
								0.5	5.0121		5.95778;
								1.0	3.61074		3.93414;
								2.0	2.36807		1.96074];
		end

		g(i) = figure();
		set(g(i), 'Position', [0 0 500 350]);
		plot(OMEGA, abs(TRANS(:,1)), 'k-', 'LineWidth', 3);
		xlabel('$He$', 'Interpreter', 'LaTeX');
		if (i == 1)
			ylabel('$\pi^+_b \slash \pi^+_a$', 'Interpreter', 'LaTeX');
			hold on;
			plot(jeonglae(:,1), jeonglae(:,2), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'b');
			h = legend('quasi-1D', 'LEE planar', 'LEE gaussian');
			set(h, 'FontSize', 14,'FontName', 'Times', 'Location', 'SouthEast');
		elseif (i == 2)
			ylabel('$\pi^+_b \slash \sigma_a$', 'Interpreter', 'LaTeX');
			hold on;
			plot(jeonglae(:,1), 2*jeonglae(:,2), 'o', 'MarkerSize', ms,'MarkerFaceColor', 'b');

			plot(jeonglae(:,1), 2*jeonglae(:,3), '>', 'MarkerSize', ms, 'MarkerFaceColor', 'r');
			h = legend('quasi-1D', 'LEE planar', 'LEE gaussian');
			set(h, 'FontSize', 14,'FontName', 'Times');
		elseif (i == 3)
			ylabel('$\pi^+_b \slash \xi_a$', 'Interpreter', 'LaTeX');
			hold on;
			plot(jeonglae(:,1), jeonglae(:,2), 'o', 'MarkerSize', ms,'MarkerFaceColor', 'b');
			plot(jeonglae(:,1), jeonglae(:,3), '>', 'MarkerSize', ms, 'MarkerFaceColor', 'r');
			h = legend('quasi-1D', 'LEE planar', 'LEE gaussian');
			set(h, 'FontSize', 14,'FontName', 'Times');
		end
		set(gca, 'FontSize', 18, 'FontName', 'Times');
		if (i == 1)
			print -depsc AcousticTransferFunction.eps
		elseif (i == 2)
			print -depsc EntropicTransferFunction.eps
		elseif (i == 3)
			print -depsc CompositionTransferFunction.eps
		end
	end%for i = 1:3
end
