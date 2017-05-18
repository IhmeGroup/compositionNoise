function[] = quickAlephCheck()
	close all;
	global mechanism;
	mechanism = 1;
	T0 = 2100;
	p0 = 1E5;
	N = 1001;
	Z = linspace(0,1,N)';
	gamma = zeros(N,1);
	for i = 1:N
		gamma(i) = returnGamma(T0,p0,Z(i));
	end
	aleph = zeros(N,1);
	dgammadZ = zeros(N,1);
	for i = 2:N-1
		dgammadZ(i) = (gamma(i+1) - gamma(i-1))/(Z(i+1) - Z(i-1));
		aleph(i) = 1/(gamma(i)*(gamma(i) - 1))*dgammadZ(i);
	end
	aleph(1) = (aleph(3) - aleph(2))/(Z(3) - Z(2))*(Z(2) - Z(1)) + aleph(2);
	aleph(N) = (aleph(N-2) - aleph(N-1))/(Z(N-2) - Z(N-1))*(Z(N) - Z(N-1)) + aleph(N-1);
	figure();
	subplot(1,2,1);
	plot(Z, gamma);
	xlabel('Z');
	ylabel('\gamma');
	set(gca, 'FontSize', 16, 'FontName', 'Times');
	subplot(1,2,2);
	plot(Z, aleph, 'LineWidth', 3);
	xlabel('Z');
	ylabel('\aleph');
	set(gca, 'FontSize', 16, 'FontName', 'Times');

	figure();
	plot(Z, dgammadZ)
	figure();
	plot(Z, 1./(gamma.*(gamma - 1)));
end
