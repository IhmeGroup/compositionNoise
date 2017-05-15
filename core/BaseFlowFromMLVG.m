function[Psibar, ubar] = BaseFlowFromMLVG(M)
	[gamma, T0, p0, Z] = returnAmbientState();
	gm1 = gamma - 1;
	gm1o2 = gm1/2;
	T = (1 + gm1o2*M*M).^(-1)*T0;
	p = (1 + gm1o2*M*M).^(-gamma/gm1)*p0;
	[Psibar] = returnPsi(T, p, Z);
	ubar = (1 + gm1o2*M*M).^(-1)*M;
end
