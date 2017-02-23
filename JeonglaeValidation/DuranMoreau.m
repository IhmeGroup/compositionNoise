function[transfer, subsol, supsol, eta, w_p, w_m, w_s, w_z, SPLINES] = DuranMoreau(M_a, M_b, M_c, Omega, runtype, suppress, SPLINY, subsol, supsol)

%	runtype is a variable to indicate what type of response is being tested. If 
%		runtype == 1 --> downstream acoustic wave
%		runtype == 3 --> entropy wave
%		runtype == 4 --> downstream acoustic wave

	addpath('../core');
	addpath('../data');

%		 1 		   2       3       4    5    6    7    8    9
%	This function drives the BVP solver for composition noise with linear velocity nozzles
%	Inputs are 
%	1) M_a = input mach number
%	2) M_b = output mach number for subsonic and transonic but unchoked flows
%	3) M_c = post-shock mach number, unused for now
%	4) Omega = dimensionless perturbation frequency
%	5) suppress = a flag to suppress generation of figures showing the spatial evolution of the system, set to true when sweeping, false when plotting a single instance/debugging
%
%	Outputs are:
%	1) transfer - a 4x2 array containing the characteristic perturbations at either end of the solution
%	2) sol - the BVP solution structure, occasionally used for debugging
%	3) eta - the spatial coordinate of the solution
%	4) w_p - (u+c) characteristic of the solution
%	5) w_m - (u-c) characteristic of the solution
%	6) w_s - entropy characteristic of the solution
%	7) w_z - composition characteristic of the solution

	global data;%This is a variable used to carry the flamelet data, without it the constant-reloading of this file slows the code down substantially
	global fuel;
	global beta;
	if (fuel <= 3) 
			data = loadFuelData(fuel); 
	end

	Omega

	epsilon = 1E-5;%This is the value of the perturbation distance (in eta) about the nozzle throat. Small values lead toward more accuracte (but more oscillatory and solower solutions)

	N1 = 401;
	N2 = 401;

	[gamma, T0, p0, Zbar] = returnAmbientState();

%	Boundary conditions imposed for the subsonic portion of the flow 
%	_p = (u+c) acoustic, _m = (u-c) acoustic, _s = entropy, _z = composition
%	Inf values are not called and intended to break the code in event of an error

	w_p_a = 0.0;
	w_m_a = Inf;
	w_s_a = 0.0;
	w_z_a = 0.0;
	w_p_b = Inf;
	w_m_b = 0.0;
	w_s_b = Inf;
	w_z_b = Inf;
	if 		(runtype == 1) w_p_a = 1;
	elseif 	(runtype == 2) error('Imposing upstream propagating acoustic wave at the inlet!');
	elseif 	(runtype == 3)	w_s_a = 1;
	elseif	(runtype == 4)	w_z_a = 1; end

%	common expressions
	gm1 = gamma - 1;%gamma - 1
	gp1 = gamma + 1;%gamma + 1
	gm1o2 = gm1/2;%(gamma - 1)/2

	if (beta == -2)
		disp('Linear Velocity Gradient Nozzle!');
		etabounds = [sqrt(gp1/2*M_a*M_a/(1+gm1o2*M_a*M_a)) sqrt(gp1/2*M_b*M_b/(1+gm1o2*M_b*M_b))];
		L = etabounds(2) - etabounds(1);
	else
		disp(strcat('Variable geometry nozzle, beta = ', num2str(beta)))
		etabounds = [-1, 1];
		L = 2;
	end


%	The param vector is used to carry thermodynamic and bc data from the driver to the ODE rhs function and BCs
	global param;
%			  1		2		3		4		5		6	7	8		9		10		11		12		13		14		15		16		17
	param = [M_a; 	M_b; 	M_c;	gamma; 	Omega; 	T0;	p0;	Zbar; 	w_p_a; 	w_m_a; 	w_s_a; 	w_z_a; 	w_p_b; 	w_m_b; 	w_s_b; 	w_z_b;  L];

	global SPLINES;
	if (~exist('SPLINY', 'var'))
		[SPLINES] = buildBaseFlowSplines();
	else
		[SPLINES] = SPLINY;
	end


%	Flags to control output
	if (suppress == true);
		plot_background			= false;
		plot_primitives 		= false;
		plot_characteristics 	= false;
		plot_invariants			= false;
	else
		close all;
		plot_background			= true;
		plot_primitives 		= true;
		plot_characteristics 	= true;
		plot_invariants			= true;
	end

	if (M_b < 1)
		subsonic = true;
		choked = false;
		shocked = false;
	else
		subsonic = false;
		if (M_c == 0)
			choked = true;
		else
			shocked = false;
		end
	end

%	Options parameter for the BVP solver
	options = [];

%	if (1 == 2)
	if (subsonic)
		disp('The flow is subsonic');
%		Compute the eta bounds (see appendix of Duran & Moreau) for the subsonic problem
		if (beta == -2)%LinVelGrad
			etabounds = [sqrt(gp1/2*M_a*M_a/(1+gm1o2*M_a*M_a)) sqrt(gp1/2*M_b*M_b/(1+gm1o2*M_b*M_b))];
		else
			etabounds = [-1,-epsilon];
		end

%		Initialize the BVP, values shouldn't matter
		if (~exist('subsol', 'var'))
			subsol 	= bvpinit(linspace(etabounds(1), etabounds(2), N1), [0.5 0.5 0.5 0.5]);
		end
		if (~exist('supsol', 'var'))
			supsol 	= bvpinit(linspace(etabounds(1), etabounds(2), N2), [0.5 0.5 0.5 0.5]);
		end

%		Solve the bvp
		subsol 	= bvp4c(@DuranMoreauODE, @SubsonicBCs, subsol, options);
%		Unpack the solution
		eta = subsol.x;		%spatial coordinate
		I_1 = subsol.y(1,:);	%mass
		I_2 = subsol.y(2,:);	%enthalpy
		I_3 = subsol.y(3,:);	%entropy
		I_4 = subsol.y(4,:);	%composition
	elseif (choked)
		disp('The flow is choked but shock free');
%		First pass (inlet -> throat)
%		options = bvpset('AbsTol', 1E-10, 'RelTol', 1E-10);
		one 	= 1.0 - epsilon;%This is the Mach number corresponding to the throat minus the perturbation distance epsilon	
%				 1		2		3		4		5		6		7		8		9		10		11		12		13		14		15		16
		param = [M_a; 	one; 	M_c;	gamma; 	Omega; 	T0;	p0;	Zbar; 	w_p_a; 	w_m_a; 	w_s_a; 	w_z_a; 	w_p_b; 	w_m_b; 	w_s_b; 	w_z_b; L];
%		Build the splines
%		[SPLINES] = buildBaseFlowSplines();
		global SPLINES;
%		Compute the corresponding spatial coordinates	
		if (beta == -2)%linVelGrad
			etabounds = [sqrt(gp1/2*M_a*M_a/(1+gm1o2*M_a*M_a)) sqrt(gp1/2*one*one/(1+gm1o2*one*one))];
		else
			etabounds = [-1,-epsilon];
		end
%		Initialize the BVP
		if (~exist('subsol', 'var'))
			subsol 	= bvpinit(linspace(etabounds(1), etabounds(2), N2), [0.5 0.5 0.5 0.5]);
		end
%		Solve the BVP
		subsol 	= bvp4c(@DuranMoreauODE, @ChokedBCs, subsol, options);
%		Unpack the BVP solution
		eta = subsol.x;		%spatial coordinate
		I_1 = subsol.y(1,:);	%mass
		I_2 = subsol.y(2,:);	%enthalpy
		I_3 = subsol.y(3,:);	%entropy
		I_4 = subsol.y(4,:);	%composition

%		Obtain the dimensionless perturbations at the throat and store them into param		
		alpha = 1/(1 + gm1o2*one*one);				% Common prefactor pre-computed for convenience
		T_t = (1 + gm1o2*one*one)^(-1)*T0;			% Temperature at the throat
		p_t = (1 + gm1o2*one*one)^(-gamma/gm1)*p0;	%pressure at the throat
		Psi_t = returnPsi(T_t, p_t, Zbar);			%Psi value at the throat

%		The perturbation invariants are related to the normalized perturbation quantities s = [p'/gamma *p0, u'/u, sigma'/c_p, xi] by the matrix P
		P = [  1	        1		 		 -1		  -Psi_t  ;
			 gm1*alpha		alpha*gm1*one*one	alpha	alpha*Psi_t;
			   0			0		  		  1			0	  ;
			   0			0		  		  0			1     ];
				
		I = [I_1(end); I_2(end); I_3(end); I_4(end)];
		s = P\I;

%		The characteristics w = (pi^+, pi^-, sigma, xi) are related to the normalized perturbation quantities s by the matrix R
		R = [1 	one 	0	0;
		   	 1 -one 	0 	0;
			 0	0		1	0;
			 0	0		0	1];

		w = R*s;

%		Re-initialize the BCs for the second pass
		w_p_a = w(1);
		w_m_a = w(2);
		w_s_a = w(3);
		w_z_a = w(4);
		w_p_b = Inf;
		w_m_b = Inf;
		w_s_b = Inf;
		w_z_b = Inf;


%		Second pass (throat -> outlet)
		one 	= 1.0 + epsilon; %Perturbation mach number on the supersonic side of the throat
%			  	 1		2		3		4		5		6		7		8		9		10		11		12		13		14		15		16 		17
		param = [one; 	M_b; 	M_c;	gamma; 	Omega; 	T0;		p0;		Zbar; 	w_p_a; 	w_m_a; 	w_s_a; 	w_z_a; 	w_p_b; 	w_m_b; 	w_s_b; 	w_z_b; 	L];
%		re-build the splines
%		[SPLINES] = buildBaseFlowSplines();
		global SPLINES;
%		Spatial bounds
		if (beta == -2)%LinVelGrad
			etabounds = [sqrt(gp1/2*one*one/(1+gm1o2*one*one)) sqrt(gp1/2*M_b*M_b/(1+gm1o2*M_b*M_b))];
		else
			etabounds = [epsilon,1];
		end 
%		Initialize
		if (~exist('supsol', 'var'))
			supsol 	= bvpinit(linspace(etabounds(1), etabounds(2), N2), [mean(I_1) mean(I_2) mean(I_3) mean(I_4)]);
		end
		supsol 	= bvp4c(@DuranMoreauODE, @SupersonicBCs, supsol, options);

%		Concatenate the subsonic and supersonic portions of the solution
		eta = [eta, supsol.x];
		I_1 = [I_1, supsol.y(1,:)];
		I_2 = [I_2, supsol.y(2,:)];
		I_3 = [I_3, supsol.y(3,:)];
		I_4 = [I_4, supsol.y(4,:)];

	elseif (shocked)
%		Shocked solution would go here...	
	end

%	Now that the full solution has been obtained, unpack it for plotting
	M_sp 		= SPLINES(1);
	Tbar_sp 	= SPLINES(2);
	pbar_sp 	= SPLINES(3);
	Psibar_sp 	= SPLINES(4);
	ubar_sp		= SPLINES(5);
	for i = 1:length(eta)
%		Compute the mean flow state from the Appendix of D&M
%		if (beta == -2)
%			etahat = eta(i).*eta(i);
%			M = sqrt((2.0/gp1)*etahat/(1 - gm1/gp1*etahat));
%		else
%			M = ppval(M_sp, eta(i));
%			Pseye = ppval(Psibar_sp, eta(i));
%		end
%		alpha = 1/(1 + gm1o2*M*M); %pre-computed for convenience
%		pee = (1 + gm1o2*M*M)^(-gamma/gm1)*p0;
%		Tee = (1 + gm1o2*M*M)^(-1)*T0;
%		Pseye = returnPsi(Tee, pee, Zbar);
		M = ppval(M_sp, eta(i));
		EM(i) = M;
		alpha = 1/(1 + gm1o2*M*M);%pre-computed for convenience
		Pseye(i) = ppval(Psibar_sp, eta(i));
		ubar(i) = ppval(ubar_sp, eta(i));
		pbar(i) = ppval(pbar_sp, eta(i));
		Tbar(i) = ppval(Tbar_sp, eta(i));
		

%		Map the invariants to the normalized perturbations s
		P = [ 	  1 		1 		  -1	  -Pseye(i);
			alpha*gm1 M*M*alpha*gm1 alpha	alpha*Pseye(i);
				  0			0 		   1	     0;
				  0			0 		   0		 1];

		I = [I_1(i); I_2(i); I_3(i); I_4(i)];
	
		s = P\I;

		phi(i) 		= s(1);
		nu(i) 		= s(2);
		sigma(i)	= s(3);
		xi(i) 		= s(4);

%		Map the normalized perturbations to the characteristics
		R = [1 	M 	0 	0;
		   	 1 -M 	0	0;
			 0	0	1	0;
			 0	0 	0 	1];
	
		w = R*s;
		w_p(i) = w(1);
		w_m(i) = w(2);
		w_s(i) = w(3);
		w_z(i) = w(4);
	end

%	Mach number vs spatial coordinate
	if (plot_background)
		figure();
		plot(eta, EM, 'b', 'LineWidth', 2);
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$M$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');
	end%(plot_primtives
	
%	Invariants vs. spatial coordinate
	if (plot_invariants)
		figure();
		plot(eta, abs(I_1), 'b', 'LineWidth', 2);
		hold on;
		plot(eta, abs(I_2), 'r', 'LineWidth', 2);
		plot(eta, abs(I_3), 'm', 'LineWidth', 2);
		plot(eta, abs(I_4), 'k', 'LineWidth', 2);
		h = legend('$I_1$', '$I_2$', '$I_3$', '$I_4$');%, '$\xi$');
		set(h, 'Interpreter','laTeX', 'FontSize', 14, 'FontName', 'Times');
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$|s_i|$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');
	end%(plot_primtives
	
%	Normalized perturbations s vs spatial coordinate
	if (plot_primitives)
		figure();
		plot(eta, abs(phi), 'b', 'LineWidth', 2);
		hold on;
		plot(eta, abs(nu), 'r', 'LineWidth', 2);
		plot(eta, abs(sigma), 'm', 'LineWidth', 2);
		plot(eta, abs(xi), 'k', 'LineWidth', 2);
		h = legend('$\phi$', '$\nu$', '$\sigma$', '$\xi$');
		set(h, 'Interpreter','laTeX', 'FontSize', 14, 'FontName', 'Times');
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$|s_i|$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');
	end%(plot_primtives
	
%	Characteristic variables vs spatial coordinate - modulus, then phase
	if (plot_characteristics)
		figure();
		plot(eta, abs(w_p), 'b', 'LineWidth', 2);
		hold on;
		plot(eta, abs(w_m), 'r', 'LineWidth', 2);
		plot(eta, abs(w_s), 'm', 'LineWidth', 2);
		plot(eta, abs(w_z), 'k', 'LineWidth', 2);
		h = legend('$w_p$', '$w_m$', '$w_s$', '$w_z$');
		set(h, 'Interpreter','laTeX', 'FontSize', 14, 'FontName', 'Times');
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$|w_i|$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		figure();
		plot(eta, atan2(imag(w_p),real(w_p)), 'b', 'LineWidth', 2);
		hold on;
		plot(eta, atan2(imag(w_m),real(w_m)), 'r', 'LineWidth', 2);
		plot(eta, atan2(imag(w_s),real(w_s)), 'm', 'LineWidth', 2);
		plot(eta, atan2(imag(w_z),real(w_z)), 'k', 'LineWidth', 2);
		h = legend('$w^p$', '$w_m$', '$w_s$','$w_z$');
		set(h, 'Interpreter','laTeX', 'FontSize', 14, 'FontName', 'Times');
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('Phase($w_i$)', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');
	end%(plot_characteristics)


	if (1 == 1)
		figure();
		plot(eta, Pseye, 'b', 'LineWidth', 2);
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$\Psi$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		figure();
		plot(eta, pbar, 'b', 'LineWidth', 2);
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$\bar{p}$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		figure();
		plot(eta, Tbar, 'b', 'LineWidth', 2);
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$\bar{T}$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		figure();
		plot(eta, ubar, 'b', 'LineWidth', 2);
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$\bar{u}$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

	end%(plot_ambient)


%	Compute the transfer matrix by storing the values at each end of the solution into the matrix transfer
	transfer = [w_p(1) w_p(end); w_m(1) w_m(end); w_s(1) w_s(end); w_z(1) w_z(end)];
end

