function[ans] = computeFluctuation(fuel)
	close all;
	if (fuel == 1)%C12H26
	data = load('../data/lowStrain/lowStrain.C12H26');
	elseif (fuel == 2)%CH4
		data = load('../data/lowStrain/lowStrain.CH4');
	elseif (fuel == 3)%H2
		data = load('../data/lowStrain/lowStrain.H2');
	end
	data = data';
	Z = data(:,1);
	T = data(:,2);
%	Z = Z_{stoich}
	Zbar = 0.0628
%	zeta is the prefactor for estimating Z'
	zeta = 1E-3;
%	Compute zprime as in the paper
	Zstd = sqrt(zeta*Zbar*(1-Zbar))
%	Based on Zbar = Zstoich and Zstd above, we can prescribe the two nonnegative parameters for the beta distribution
	alpha = -(Zbar*(Zbar^2 - Zbar + Zstd*Zstd))/(Zstd*Zstd)
	beta = (Zbar^3 - 2*Zbar^2 + Zbar*Zstd^2 + Zbar - Zstd^2)/Zstd^2
	Tbar = 2100;
%	Beta PDF normalization factor. Alpha and beta are both large, so this is difficult to evaluate numerically. Instead build the un-normalized pdf, integrate it and normalize. (I'm a numericist!)
%	denom = gamma(alpha)*gamma(beta)/gamma(alpha + beta)
	N = length(T)
	for i = 1:N
		b(i) = Z(i).^(alpha - 1).*(1 - Z(i)).^(beta - 1);
	end
	b = b/trapz(Z, b);
%	Plot the beta PDF to be sure it looks like science
	figure();
	plot(Z, b, 'LineWidth', 2);
	xlabel('Z');
	ylabel('\beta(Z)');
%	Compute the integrand? in expression 36. 
	for i = 1:N
%		pdf(i) = (T(i) - Tbar).^2.*Z(i).^(alpha - 1).*(1-Z(i)).^(beta - 1)/denom;
		pdf(i) = (T(i) - Tbar).^2.*b(i);
	end
%	Numerically integrate it and take the square root
	Tprime = sqrt(trapz(Z, pdf))
%	Output the ratio	
	ans = Zstd/Tprime
%	ans = Tprime;

end
