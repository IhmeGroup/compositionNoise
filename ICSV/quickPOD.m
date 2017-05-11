function[s_fit, z_fit] = quickPOD(pi,sigma,xi)
	close all;
	data = load('processedFlucs.mat');
	v = data.v_r;
	p = data.p_r;
	rho = data.rho_r;
	z = data.Z_r;
	pbar = data.pbar;
	Tbar = data.Tbar;
	rhobar = data.rhobar;

%	Nozzle entrance condition
	gamma = 1.4;
	Tbar_1 = mean(Tbar)
	pbar_1 = mean(pbar)
	rhobar_1 = mean(rhobar);
	Psibar_1 = -21;
	M_1 = 0.01;
	M_2 = 1.5;
	cbar_1 = sqrt(1.4*287*Tbar_1);
	c_p = 1004;

	[Nsnaps, Npts] = size(p)

	pi = p./(gamma*pbar_1) + v./cbar_1;
	sigma = p./(gamma*pbar_1) - rho/rhobar_1 - Psibar_1*z;
	xi = z;

	pibar = mean(pi);
	sigmabar=  mean(sigma);
	xibar = mean(xi);

	for i = 1:Nsnaps
		pi(i,:) = pi(i,:) - pibar;
		sigma(i,:) = sigma(i,:) - sigmabar;
		xi(i,:) = xi(i,:) - xibar;
	end

	size(pi)
	size(mean(pi))


	[phi_p, lmd_p, r_p] = pod(pi');
	[phi_s, lmd_s, r_s] = pod(sigma');
	[phi_z, lmd_z, r_z] = pod(xi');

	h1 = figure();
	plot(phi_p(:,1:5), 'LineWidth', 2);
	legend('1', '2', '3','4', '5');
	title('Pi');

	h2 = figure();
	plot(phi_s(:,1:5), 'LineWidth', 2);
	legend('1', '2', '3','4', '5');
	title('Sigma');

	h3 = figure();
	plot(phi_z(:,1:5), 'LineWidth', 2);
	legend('1', '2', '3','4', '5');
	title('Xi');

	r = linspace(0,1,Npts);
	size(r)
	size(phi_s(:,1))

	s_fit = polyfit(r', phi_s(:,1), 6);
	z_fit = polyfit(r', phi_z(:,1), 6);

	R = linspace(0,1,1001)';
	S = zeros(1001,1);
	Z = zeros(1001,1);
	size(S)
	size(Z)
	size(R)
	for i = 1:1001
		S(i) = polyval(s_fit, R(i));
		Z(i) = polyval(z_fit, R(i));
	end

	figure(h2);
	hold on;
	plot(51*R, S, 'r--', 'LineWidth', 3);
	figure(h3);
	hold on;
	plot(51*R, Z, 'r--', 'LineWidth', 3);


end
