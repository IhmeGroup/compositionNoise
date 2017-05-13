
	elseif ((beta <= 1) && (beta >= -1))
		etabounds = [-1,1];
		eta = linspace(etabounds(1), etabounds(2), N);
		AoverAstar = zeros(N,1);
		A1 = zeros(N,1);
		A2 = zeros(N,1);
		A3 = zeros(N,1);
		for i = 1:N
			if (eta(i) < 0)
				A1(i) = 1 - eta(i);
				A3(i) = 1 + sqrt(1 - (eta(i) + 1).^2);
			else
				A1(i) = 1 + eta(i);
				A3(i) = 1 + sqrt(1 - (eta(i) - 1).^2);
			end
			A2(i) = 2 - sqrt(1 - eta(i).^2);
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
	end%Geometry study instead of linear vel nozzle	
