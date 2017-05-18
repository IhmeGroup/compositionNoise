function[] = testBaseFlow(M_a, M_b, M_c)
	close all;
	SPLINES = [M_a, M_b, M_c];
	N = 1001;
	eta = linspace(0,1,N)';
	M = zeros(N,1);
	for i = 1:N
		M(i) = MFromEtaDMSC(eta(i), [M_a M_b M_c]);
	end
	figure();
	plot(eta, M, 'LIneWidth', 3);
end
