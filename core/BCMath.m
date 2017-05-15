function[soln] = BCMath()
	Mb = sym('Mb', 'positive');
	Mc = sym('Mc', 'positive');

	pbp = sym('pbp');
	pbm = sym('pbm');
	sb = sym('sb');
	pcp = sym('pcp');
	pcm = sym('pcm');
	sc = sym('sc');
	gamma = sym('gamma', 'real');
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;

	alphab = 1/(1 + gm1o2*Mb*Mb);
	alphac = 1/(1 + gm1o2*Mc*Mc);

	ConsMassLH 	= 1/2*(1/Mb+1)*pbp - 1/2*(1/Mb-1)*pbm - sb;
	ConsMassRH 	= 1/2*(1/Mc+1)*pcp - 1/2*(1/Mc-1)*pcm - sc;
	Eq1 = ConsMassLH - ConsMassRH

	ConsMomLH	= alphab*gm1*(Mb+1)/2*pbp - alphab*gm1*(Mb-1)/2*pbm + alphab*sb;
	ConsMomRh 	= alphac*gm1*(Mc+1)/2*pcp - alphac*gm1*(Mc-1)/2*pcm + alphac*sc;
	Eq2	= ConsMomLH - ConsMassRH

	soln = simplify(solve(Eq1, pcp));
	Eq2 = subs(Eq2, pcp, soln);
	Eq2 = simplify(Eq2)
	soln = solve(Eq2, sc)
%	soln = solve(Eq1, Eq2, pcp, sc);

	Mb = 1.5;
	pbm = 0.9227
	pbp = 1.1680;
	sb	= 0;

	gamma = returnAmbientState();
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;
	alphab = 1/(1 + gm1o2*Mb*Mb);

	ConsMassLH 	= 1/2*(1/Mb+1)*pbp - 1/2*(1/Mb-1)*pbm - sb
	ConsMomLH	= alphab*gm1*(Mb+1)/2*pbp - alphab*gm1*(Mb-1)/2*pbm + alphab*sb



end
