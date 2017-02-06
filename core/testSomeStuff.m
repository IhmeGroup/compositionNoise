function[] = testSomeStuff(beta)
	close all;
	N = 101;
	[gamma, T0, p0, Zbar] = returnAmbientState()
		etabounds = [-1,1];
		eta = linspace(etabounds(1), etabounds(2), N)';
		AoverAstar = zeros(N,1);
		A1 = zeros(N,1);
		A2 = zeros(N,1);
		A3 = zeros(N,1);
		M 		= zeros(N,1);
		Tbar 	= zeros(N,1);
		pbar 	= zeros(N,1);
		Psibar	= zeros(N,1);
		ubar 	= zeros(N,1);
		gm1 = gamma - 1;
		gp1 = gamma + 1;
		gm1o2 = gm1/2;
		gp1o2 = gp1/2;
		for i = 1:N
			if (eta(i) < 0)
				A1(i) = 1 - eta(i);
				A2(i) = 1 + sqrt(1 - (eta(i) + 1).^2);
			else
				A1(i) = 1 + eta(i);
				A2(i) = 1 + sqrt(1 - (eta(i) - 1).^2);
			end
			A3(i) = 2 - sqrt(1 - eta(i).^2);
		end
		if (beta < 0) 
			delta = -beta;	
			AoverAstar = (1-delta)*A1 + delta*A2;
		else
			AoverAstar = (1-beta)*A1 + beta*A3;
		end
		for i = 1:N
			if (eta(i) < 0) 
				Mguess = 0.9;
			else
				Mguess = 1.05;
			end
			options = optimset('Display', 'off');
			M(i) = fsolve(@MfromA, Mguess, options, AoverAstar(i));
			M2 = M(i)*M(i);	
			Tbar(i) 	= (1 + gm1o2*M2)^(-1)*T0;
			pbar(i) 	= (1 + gm1o2*M2)^(-gamma/gm1)*p0;
			Psibar(i) 	= returnPsi(Tbar(i), pbar(i), Zbar);
			ubar(i) 	= M(i)*sqrt(Tbar(i)/T0);
		end
	plot(eta, A1, 'b-', 'LineWidth', 4);
	hold on;
	plot(eta, A2, 'r--', 'LineWidth', 4);
	plot(eta, A3, 'g:', 'LineWidth', 4);
	plot(eta, AoverAstar, 'k-', 'LineWidth', 4);
	ylim([0, 2]);
	legend('\beta = 0', '\beta = 1', '\beta = -1');

	figure();
	subplot(2,3,1);
	plot(eta, M)
	title('Mach');
	subplot(2,3,2);
	plot(eta, M2)
	title('M2');
	subplot(2,3,3);
	size(eta)
	size(Tbar)
	plot(eta, Tbar)
	title('T');
	subplot(2,3,4);
	plot(eta, pbar)
	title('p');
	subplot(2,3,5);
	plot(eta, Psibar)
	title('Psi');
	subplot(2,3,6);
	plot(eta, ubar)
	title('u');
end%function
