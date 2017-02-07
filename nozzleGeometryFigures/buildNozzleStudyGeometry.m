function[] = buildNozzleStudyGeometry()
%	This is a function that builds the spatial coordinate vs. cross-sectional area and spatial coordinate vs Mach number figures for the geometry study at the values beta = -1, 0, and 1

	close all;
	lw = 4;
	fs = 18;


%	Compute the cross sectional areas	
	N = 500;
	x = linspace(-1,1,N)';
	A1 = [linspace(2,1,N/2)'; linspace(1,2,N/2)'];
	A2 = 2-sqrt(1 - (x).^2);
	A3 = sqrt(1 - (x).^2);
	A3 = [A3(N/2+1:end); A3(1:N/2)] + 1;
%	shift x to keep the throat at 1
	x = x + 1;

%	Rescale x so L =1 for all nozzles
	ex = x - min(x);
	ex = ex./max(ex);

%	plot the cross sectional area vs. eta
	h = figure();
	set(h, 'Position', [0 0 800 300]);
	subplot('Position', [0.08 0.15 0.4 0.75]);
	plot(ex, A1, 'k-', 'LineWidth', lw);
	hold on;
	plot(ex, A2, 'b--','LineWidth', lw);
	plot(ex, A3, 'r:', 'LineWidth', lw);
	plot([.5 .5], [0, 1], 'k--', 'LineWidth', lw);

%	and format it
	ylim([-0,2]);
	xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
	ylabel('$A(\eta) \slash A^*$', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
	set(gca, 'FontSize', fs, 'FontName', 'Times');
	title('$(a)$', 'FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');


%	Compute the Mach number using an iterative solver
	M1 = zeros(N,1);
	M2 = zeros(N,1);
	M3 = zeros(N,1);

	for i= 1:N
		disp(strcat(num2str(round(1000*i/N)/10), '% complete'));
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

%	Rescale x so L =1 for all nozzles
	x = x - min(x);
	x = x./max(x);
		
%	Plot the profile
	subplot('Position', [0.58 0.15 0.4 0.75]);
	plot(x, M1, 'k-', 'LineWidth', lw);
	hold on;
	plot(x, M2, 'b--', 'LineWidth', lw);
	plot(x, M3, 'r:', 'LineWidth', lw);
	xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
	ylabel('$M(\eta)$', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
%	h = legend('$\beta = 0$', '$\beta = 1$', '$\beta = -1$','Location', 'Southeast');
%	set(h, 'FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	title('$(b)$', 'FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	set(gca, 'FontSize', fs, 'FontName', 'Times');

	print -depsc NozzleGeometry.eps
end%buildNozzleGeometry
