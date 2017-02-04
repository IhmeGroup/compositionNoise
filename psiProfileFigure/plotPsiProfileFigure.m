function[] = plotPsiProfileFigure()
	close all;

	addpath('../core');
	addpath('../data');

	global data;%This is a variable used to carry the flamelet data, without it the constant-reloading of this file slows the code down substantially
	global fuel, param, SPLINES, beta;
	beta = -2;%Use linear vel grad nozzle
	data = loadFuelData(fuel);

	[gamma, T0, p0, Zbar] = returnAmbientState();

%	Boundary conditions imposed for the subsonic portion of the flow 
%	_p = (u+c) acoustic, _m = (u-c) acoustic, _s = entropy, _z = composition
%	Inf values are not called and intended to break the code in event of an error
	w_p_a = Inf;
	w_m_a = Inf;
	w_s_a = Inf;
	w_z_a = Inf;
	w_p_b = Inf;
	w_m_b = Inf;
	w_s_b = Inf;
	w_z_b = Inf;
	M_c = Inf;
	Omega = 0;
	M_a = 0.29;
	
	lw = 4;

%	common expressions
	gm1 = gamma - 1;%gamma - 1
	gp1 = gamma + 1;%gamma + 1
	gm1o2 = gm1/2;%(gamma - 1)/2

	N = 401;
	ETA = zeros(N, 3);
	PSI = zeros(N, 3);
	for run = 1:3
		if (run == 1) M_b = 0.89;
		elseif (run == 2) M_b = 1.02;
		elseif (run == 3) M_b = 1.50; end

		etabounds = [sqrt(gp1/2*M_a*M_a/(1+gm1o2*M_a*M_a)) sqrt(gp1/2*M_b*M_b/(1+gm1o2*M_b*M_b))];
		L = etabounds(2) - etabounds(1);

%		The param vector is used to carry thermodynamic and bc data from the driver to the ODE rhs function and BCs
%				  1		2		3		4		5		6	7	8		9		10		11		12		13		14		15		16
		param = [M_a; 	M_b; 	M_c;	gamma; 	Omega; 	T0;	p0;	Zbar; 	w_p_a; 	w_m_a; 	w_s_a; 	w_z_a; 	w_p_b; 	w_m_b; 	w_s_b; 	w_z_b; L];

		[SPLINES] = buildBaseFlowSplines();

		Psibar_sp 	= SPLINES(4);
	
		eta = linspace(etabounds(1), etabounds(2), N)';
		for i = 1:N
			psi(i) = ppval(Psibar_sp, eta(i));
		end

		eta = eta - min(eta);
		eta = eta./max(eta);

		ETA(:,run) = eta;
		PSI(:,run) = psi;
	end

	figure();
	plot(ETA(:,1), PSI(:,1), 'b-', 'LineWidth', lw);
	hold on;
	plot(ETA(:,1), PSI(:,2), 'r--', 'LineWidth', lw);
	plot(ETA(:,1), PSI(:,3), 'k:', 'LineWidth', lw);
	xlabel('$\eta$', 'Interpreter', 'laTeX', 'FontSize', 24, 'FontName', 'Times');
	ylabel('$\Psi(\eta)$', 'FontSize', 24, 'FontName', 'Times', 'Interpreter','LaTeX');
	set(gca, 'FontSize', 24, 'FontName', 'Times');
	grid on;

	print -depsc PsiProfile.eps
end%function
