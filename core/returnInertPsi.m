%		 1 	  2   3     4     5     6						   1  2  3     4     5     6     7     8     9  10 11  12    13	
function[psi, cp, hctr, sctr, gctr, gammactr] = returnInertPsi(T, p, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW, Rbar, cpfixed)
	[~, ~, ~, glft] = returnSpeciesProperties(T, p, Ylft, a, A, MW);
	[cp, hctr, sctr, gctr]	= returnSpeciesProperties(T, p, Yctr, a, A, MW);
	[~, ~, ~, grgt] = returnSpeciesProperties(T, p, Yrgt, a, A, MW);

	cv = cp - Rbar;

	gammactr = cp/cv;

	dZ1 = Zctr - Zlft;
	dZ2 = Zctr - Zrgt;

	dg1 = gctr - glft;
	dg2 = gctr - grgt;

	denom = (dZ1*(-dZ2^2 + dZ1*dZ2));
	
	dGdZ = -(-dg2*dZ1^2 + dg1*dZ2^2)/denom;
%	dGdZ = (grgt - glft)/(Zrgt - Zlft);
	psi = 1/(cpfixed*T)*dGdZ;
end
