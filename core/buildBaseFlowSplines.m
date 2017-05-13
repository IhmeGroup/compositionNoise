function[SPLINES] = buildBaseFlowSplines(M_a, M_b, N)
	[gamma, T0, p0, Zbar] = returnAmbientState();

%	Precompute for speed I guess	
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;

%	Pre-allocate Memory
	M 		= zeros(N,1);
	Tbar 	= zeros(N,1);
	pbar 	= zeros(N,1);
	Psibar	= zeros(N,1);
	ubar 	= zeros(N,1);

	etabounds = [sqrt(gp1/2*M_a*M_a/(1+gm1o2*M_a*M_a)) sqrt(gp1/2*M_b*M_b/(1+gm1o2*M_b*M_b))];
	eta = linspace(etabounds(1), etabounds(2), N);
	for i = 1:N
		etahat = eta(i)*eta(i);
		M(i) = sqrt((2/gp1)*etahat/(1-gm1/gp1*etahat));
		M2 = M(i)*M(i);	
		Tbar(i) 	= (1 + gm1o2*M2)^(-1)*T0;
		pbar(i) 	= (1 + gm1o2*M2)^(-gamma/gm1)*p0;
		Psibar(i) 	= returnPsi(Tbar(i), pbar(i), Zbar);
		ubar(i) 	= M(i)*sqrt(Tbar(i)/T0);
	end
	eta = (eta - min(eta))./(max(eta) - min(eta));

	M_sp 		= spline(eta,M);
	Tbar_sp 	= spline(eta, Tbar);
	pbar_sp		= spline(eta, pbar);
	Psibar_sp	= spline(eta, Psibar);
	ubar_sp		= spline(eta, ubar);

	global SPLINES;
%				1		2		3		4			5	
	SPLINES	= [M_sp; Tbar_sp; pbar_sp; Psibar_sp; ubar_sp];
	size(SPLINES);
	disp('build SPLINES');
end
