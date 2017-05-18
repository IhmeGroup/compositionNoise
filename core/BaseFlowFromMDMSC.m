function[Psibar, ubar] = BaseFlowFromMDMSC(eta, M)
	[gamma, T0, p0, Z] = returnAmbientState();
	gm1 = gamma - 1;
	gm1o2 = gm1/2;
	T = (1 + gm1o2*M*M).^(-1)*T0;
	if (eta < 0.5008)
		p = (1 + gm1o2*M*M).^(-gamma/gm1)*p0;
	else
		p = (1 + gm1o2*M*M).^(-gamma/gm1)*p0*0.919;
	end
	[Psibar] = returnPsi(T, p, Z);
	ubar = (1 + gm1o2*M*M).^(-1)*M;
end
