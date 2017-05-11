function[] = generateRawDataForJeonglae()
	close all;
	N = 101;
	addpath('../core');
	global fuel;
	fuel = 5;
	[gamma, T0, p0, Zbar] = returnAmbientState();
	p0
	gm1 = gamma - 1;
	gm1o2 = gm1/2;
	M = linspace(0,2,N)';
	u = zeros(N,1);
	p = zeros(N,1);
	T = zeros(N,1);
	pstar = zeros(N,1);
	Tstar = zeros(N,1);
	Psi = zeros(N,1);
	for i = 1:N
		T(i) = T0*(1 + gm1o2*M(i)^2).^(-1);
		Tstar(i) = (1 + gm1o2*M(i)^2).^(-1);
		p(i) = p0*(1 + gm1o2*M(i).^2).^(-gamma/gm1);
		pstar(i) = (1 + gm1o2*M(i).^2).^(-gamma/gm1);
		[Psi(i),~,Ratio(i)] = returnPsi(T(i), p(i), Zbar);
	end%for i = 1:N

	figure();
	plot(M, Ratio);
	legend('RATIO');
	xlabel('M');
	set(gca, 'FontSize', 16);

	figure();
	plot(M, [Tstar, pstar], 'LineWidth', 3);
	legend('T/T0', 'p/p0');
	xlabel('M');
	set(gca, 'FontSize', 16);

	figure();
	plot(M, T, 'LineWidth', 3);
	xlabel('M');
	ylabel('T [K]');
	set(gca, 'FontSize', 16);

	figure();
	plot(M, p./1E5, 'LineWidth', 3);
	xlabel('M');
	ylabel('p [bar]');
	set(gca, 'FontSize', 16);

	figure();
	plot(M, Psi, 'LineWidth', 3);
	xlabel('M');
	ylabel('\Psi [-]');
	set(gca, 'FontSize', 16);
	
	save 'MeanFlow.mat' M u p T pstar Tstar Psi;


end%function
