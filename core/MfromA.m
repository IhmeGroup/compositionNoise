function [res] = MfromA(M,A)
%Returns the residual of the mach-normalized area relation for help inverting the equation to obtain M from A
	gamma = 1.4;
	gm1 = gamma - 1;
	gp1o2 = (gamma+1)/2;
	gm1o2 = (gamma-1)/2;
	res = gp1o2.^(-gp1o2/gm1)*(1+gm1o2*M*M).^(gp1o2/gm1)/M - A;
end
