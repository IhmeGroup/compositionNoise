function[] = buildNozzleStudyGeometry()
%	This is a function that builds the spatial coordinate vs. cross-sectional area and spatial coordinate vs Mach number figures for the geometry study at the values beta = -1, 0, and 1

	close all;


%	Compute the cross sectional areas	
	N = 1000;
	x = linspace(-1,1,N)';
	A1 = [linspace(2,1,N/2)'; linspace(1,2,N/2)'];
	A2 = 2-sqrt(1 - (x).^2);
	A3 = sqrt(1 - (x).^2);
	A3 = [A3(N/2+1:end); A3(1:N/2)] + 1;
%	shift x to keep the throat at 1
	x = x + 1;

%	plot the cross sectional area vs. eta
	figure();
	plot(x, A1, 'k-', 'LineWidth', 2);
	hold on;
	plot(x, A2, 'b--','LineWidth', 2);
	plot(x, A3, 'r:', 'LineWidth', 2);
	plot([1 1], [0, 1], 'k--', 'LineWidth', 2);

%	and format it
	ylim([-0,2]);
	xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$A(\eta) \slash A^*$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName', 'Times');
	h = legend('$\beta = 0$', '$\beta = 1$', '$\beta = -1$','Location', 'Southeast');
	set(h, 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'LaTeX');


%	Compute the Mach number using an iterative solver
	M1 = zeros(N,1);
	M2 = zeros(N,1);
	M3 = zeros(N,1);

	for i= 1:N
%		Seed the iterative solver with a subsonic guess if x < 1 or a supersonic guess if x > 1
		if (x(i) < 1) Mguess = 0.5;
		else Mguess = 1.5;
		end

		if (i == N/2)
			M1(i) = 1;
			M2(i) = 1;
			M3(i) = 1;
		else
			M1(i) = fsolve(@MfromA, Mguess, [], A1(i));
			M2(i) = fsolve(@MfromA, Mguess, [], A2(i));
			M3(i) = fsolve(@MfromA, Mguess, [], A3(i));
		end
	end
		
%	Plot the profile
	figure();
	plot(x, M1, 'k-', 'LineWidth', 2);
	hold on;
	plot(x, M2, 'b--', 'LineWidth', 2);
	plot(x, M3, 'r:', 'LineWidth', 2);
	xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$\bar{M}(\eta)$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	h = legend('$\beta = 0$', '$\beta = 1$', '$\beta = -1$','Location', 'Southeast');
	set(h, 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'LaTeX');

end%buildNozzleGeometry
