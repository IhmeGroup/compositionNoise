%			1		 2		 3		 4		  5	           		   1    2    3     4      5          6         7        8        9     10       11      12
function[transfer, subsol, supsol, shksol, SPLINES] = DuranMoreau(M_a, M_b, M_c, Omega, runtype, mechansm, suppress, flaggoe, SPLINES, subsol, supsol, shksol)
%	This function drives the BVP solver for composition noise with linear velocity nozzles
%	Inputs are 
%	1) M_a = input mach number
%	2) M_b = output mach number for subsonic and transonic but unchoked flows
%	3) M_c = post-shock mach number, unused for now
%	4) Omega = dimensionless perturbation frequency
%	5) runtype is a variable to indicate what type of response is being tested. If 
%		runtype == 1 --> downstream acoustic wave
%		runtype == 3 --> entropy wave
%		runtype == 4 --> downstream acoustic wave
%	6) mechanism identifies which fuel/application is being considered here
%		mechanism == 1 --> C12H26 (dodecane)
%		mechanism == 2 --> CH4	(methane)
%		mechanism == 3 --> H2 (Hydrogen)
%		mechanism == 4 --> Inerts
%		mechanism == 5 --> MW
%	7) suppress = a flag to suppress generation of figures showing the spatial evolution of the system, set to true when sweeping, false when plotting a single instance/debugging
%	8) flaggo is a way to switch between a splined background flow and a locally computed background flow
%		flaggo == 0 --> standard spline (good for shock-free applications with linear velocity gradient geometry)
%		flaggo == 1 --> Locally computed M(A) (i.e. no Gibbs) for linear velocity gradient nozzle. Serves as check on 
%						flaggo == 0 (though much slower) and also good for generic shocked applications
%		flaggo == 2 --> Bakke's transonic nozzle (need local M(A) solution b/c of shock and Gibbs phenomenon going batshit here)
%		flaggo == 3 --> Smoothly splined, but with the stupid convex to concave circle geometry to assess geometric sensitivity
%	9)  SPLINES contains the spline object to allow for recycling in parameter sweep applications
%	10) subsol contains the subsonic solution object to allow for recyling in parameter sweep applications
%	12) supsol contains the supersonic solution object to allow for recyling in parameter sweep applications
%	12) shksol contains the post-shock solution object to allow for recyling in parameter sweep applications
%
%	Outputs are:
%	1) transfer - a 4x2 array containing the characteristic perturbations at either end of the solution
%	2) subsol - the subsonic portion of the BVP solution structure, used for solution recyling
%	3) supsol - the supersonic portion of the BVP solution structure, used for solution recyling
%	4) shksol - the post-shock portion of the BVP solution structure, used for solution recyling
%	5) SPLINES - the baseflow structure for solution recyling

%	Global variables are cludgy as $#@!, but this is simpler than making the whole code object oriented
	global data;%This is a variable used to carry the flamelet data, without it the constant-reloading of this file slows the code down substantially
	global mechanism;
	mechanism = mechansm;
	global flaggo;
	flaggo = flaggoe;
	global epsilon;

%	Adda a bunch of directories to the path so we can find crap
	addpath('../core');
	addpath('../data');
	addpath('../speciesProps');

	Omega
%	Solver parameters - hard-coding some stuff so that the inputs don't get (more) verbose
	epsilon = 1E-4;%This is the value of the perturbation distance (in eta) about the nozzle throat. Small values lead toward more accuracte (but more oscillatory and solower solutions)
	N1 = 401;
	N2 = 401;
	N3 = 401;
	Nsplines = 101;

%	Preload the flamelet data for Z gradients
	loadSpeciesData();

%	Grab the stagnation condition for this run
	[gamma, T0, p0, Zbar] = returnAmbientState();

%	Precomputing some common expressions
	gm1 = gamma - 1;%gamma - 1
	gp1 = gamma + 1;%gamma + 1
	gm1o2 = gm1/2;%(gamma - 1)/2
	gp1o2 = gp1/2;%(gamma - 1)/2

%	Flags to control output
	if (suppress == true);
		plot_background			= false;
		plot_primitives 		= false;
		plot_characteristics 	= false;
		plot_invariants			= false;
		plot_psi				= false;
	else
		close all;
		plot_background			= true;
		plot_primitives 		= true;
		plot_characteristics 	= true;
		plot_invariants			= true;
		plot_psi				= true;
	end

%	Determine the flow regime
	choked 		= false;
	shocked 	= false;
	subsonic 	= false;
	supersonic 	= false;
	diffusor 	= false;
	if (M_a < 1) && (M_b < 1)
		subsonic 	= true;
		disp('The flow is subsonic throughout the entire domain');
	elseif (M_a < 1) && (M_b > 1) && (M_c == 0)
		choked 		= true;
		disp('The flow is choked but shock free');
	elseif (M_a < 1) && (M_b > 1) && (M_c ~= 0)
		shocked 	= true;
		disp('The flow is choked and experiences a downstream shock');
	elseif (M_a > 1) && (M_b > 1)
		supersonic 	= true;
		disp('The flow is completely supersonic');
	else
		error('Unsupported configuration!');
	end

%	Options parameter for the BVP solver
	options = bvpset('AbsTol', 1E-8, 'RelTol', 1E-8);
%	options = bvpset('AbsTol', 1E-5, 'RelTol', 1E-5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		SUBSONIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	if (subsonic)
%		Do something about BCs		
		if 		(runtype == 1) w_l = [1; NaN; 0; 0]; w_r = [NaN; 0; NaN; NaN];
		elseif 	(runtype == 2) error('Imposing upstream propagating acoustic wave at the inlet!');
		elseif 	(runtype == 3) w_l = [0; NaN; 1; 0]; w_r = [NaN; 0; NaN; NaN];	
		elseif	(runtype == 4) w_l = [0; NaN; 0; 1]; w_r = [NaN; 0; NaN; NaN];
		end

		if (flaggo == 0)
			if (~exist('SPLINES'))
%				Fresh start for splined linear velocity gradient				
				disp('Splined base flow computation');
				SPLINES = buildBaseFlowSplines(M_a, M_b, Nsplines);
			else
				disp('Recycled splined base flow computation');
			end
		elseif (flaggo == 1)
%				Locally based computation - this is an ugly variable container					
				disp('Local base flow computation');
				SPLINES = [M_a, M_b, M_c];
		else
			error('Shoulda written more code');
		end

		ODECaller = @(eta, I) DuranMoreauODE(eta, I, Omega, flaggo, SPLINES);
		BCCaller = @(I_a, I_r) SubsonicBCs(I_a, I_r, w_l, w_r, 0, 1, flaggo, SPLINES);

%		Initialize the BVP, values shouldn't matter
		if (~exist('subsol', 'var'))
			subsol 	= bvpinit(linspace(0, 1, N1), [0.5 0.5 0.5 0.5]);
		end

%		Solve the bvp
		subsol 	= bvp4c(ODECaller, BCCaller, subsol, options);
		supsol 	= NaN;
		shksol 	= NaN;

%		Unpack the solution
		eta = subsol.x;		%spatial coordinate
		I_1 = subsol.y(1,:);	%mass
		I_2 = subsol.y(2,:);	%enthalpy
		I_3 = subsol.y(3,:);	%entropy
		I_4 = subsol.y(4,:);	%composition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		SUPERSONIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	elseif (supersonic)
%		Do something about BCs		
		if 		(runtype == 1) w_l = [1; 0; 0; 0]; w_r = [NaN; NaN; NaN; NaN];
		elseif 	(runtype == 2) w_l = [0; 1; 0; 0]; w_r = [NaN; NaN; NaN; NaN];
		elseif 	(runtype == 3) w_l = [0; 0; 1; 0]; w_r = [NaN; NaN; NaN; NaN];	
		elseif	(runtype == 4) w_l = [0; 0; 0; 1]; w_r = [NaN; NaN; NaN; NaN];
		end

%		Compute the eta bounds (see appendix of Duran & Moreau) for the supersonic problem
		if (~exist('SPLINES'))
			if (flaggo == 0) 
%				Fresh start for splined linear velocity gradient				
				disp('Splined base flow computation');
				SPLINES = buildBaseFlowSplines(M_a, M_b, Nsplines);
			else
%				Locally based computation - this is an ugly variable container					
				disp('Local base flow computation');
				SPLINES = [M_a, M_b, M_c];
			end
		end

		ODECaller = @(eta, I) DuranMoreauODE(eta, I, Omega, flaggo, SPLINES);
		BCCaller = @(I_a, I_r) SupersonicBCs(I_a, I_r, w_l, w_r, 0, 1, flaggo, SPLINES);

%		Initialize the BVP, values shouldn't matter
		if (~exist('supsol', 'var'))
			supsol 	= bvpinit(linspace(0, 1, N2), [0.5 0.5 0.5 0.5]);
		end

%		Solve the bvp
		subsol	= NaN;
		supsol 	= bvp4c(ODECaller, BCCaller, supsol, options);
		shksol	= NaN;

%		Unpack the solution
		eta = supsol.x;		%spatial coordinate
		I_1 = supsol.y(1,:);	%mass
		I_2 = supsol.y(2,:);	%enthalpy
		I_3 = supsol.y(3,:);	%entropy
		I_4 = supsol.y(4,:);	%Mixture Fraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		CHOKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	elseif (choked)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%			First Pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		First pass (inlet -> throat)
%		Do something about BCs		
		if 		(runtype == 1) w_l = [1; NaN; 0; 0]; w_r = [NaN; NaN; NaN; NaN];
		elseif 	(runtype == 2) error('Imposing upstream propagating acoustic wave at the inlet!');
		elseif 	(runtype == 3) w_l = [0; NaN; 1; 0]; w_r = [NaN; NaN; NaN; NaN];	
		elseif	(runtype == 4) w_l = [0; NaN; 0; 1]; w_r = [NaN; NaN; NaN; NaN];
		end

%		Compute the eta bounds (see appendix of Duran & Moreau) for the supersonic problem
		if (~exist('SPLINES'))
			if (flaggo == 0) 
%				Fresh start for splined linear velocity gradient				
				disp('Splined base flow computation');
				SPLINES = buildBaseFlowSplines(M_a, M_b, Nsplines);
			elseif (flaggo == 1)
%				Locally based computation - this is an ugly variable container					
				disp('Local base flow computation');
				SPLINES = [M_a, M_b, M_c];
			else
				error('Need to address other flaggo cases here');
			end
		end

		x_a 	= sqrt(gp1/2*M_a*M_a/(1+gm1o2*M_a*M_a));
		x_b 	= sqrt(gp1/2*M_b*M_b/(1+gm1o2*M_b*M_b));
		x_1 	= sqrt(gp1/2*1*1/(1+gm1o2*1*1));
		x_1		= (x_1 - x_a)./(x_b - x_a);
		etabounds = [0, x_1-epsilon];
		ODECaller = @(eta, I) DuranMoreauODE(eta, I, Omega, flaggo, SPLINES);
		BCCaller = @(I_a, I_r) ChokedBCs(I_a, I_r, w_l, w_r, etabounds(1), etabounds(2), flaggo, SPLINES);

%		Initialize the BVP
		if (~exist('subsol', 'var'))
			subsol 	= bvpinit(linspace(etabounds(1), etabounds(2), N2), [0.5 0.5 0.5 0.5]);
		end
%		Solve the BVP
		subsol 	= bvp4c(ODECaller, BCCaller, subsol, options);

%		Unpack the BVP solution
		eta = subsol.x;		%spatial coordinate
		I_1 = subsol.y(1,:);	%mass
		I_2 = subsol.y(2,:);	%enthalpy
		I_3 = subsol.y(3,:);	%entropy
		I_4 = subsol.y(4,:);	%composition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%			Second Pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		Do something about the BCs
		I_end 	= subsol.y(:,end);
		[w_l] = charsFromInvrnts(eta(end), I_end, flaggo, SPLINES);
		w_r		= [NaN; NaN; NaN; NaN];

%		Compute the eta bounds (see appendix of Duran & Moreau) for the supersonic problem
		etabounds = [x_1+epsilon, 1];
		ODECaller = @(eta, I) DuranMoreauODE(eta, I, Omega, flaggo, SPLINES);
		BCCaller = @(I_a, I_r) SupersonicBCs(I_a, I_r, w_l, w_r, etabounds(1), etabounds(2), flaggo, SPLINES);

%		Initialize
		if (~exist('supsol', 'var'))
			supsol 	= bvpinit(linspace(etabounds(1), etabounds(2), N2), [mean(I_1) mean(I_2) mean(I_3) mean(I_4)]);
		end
		supsol 	= bvp4c(ODECaller, BCCaller, supsol, options);
		shksol	= NaN;

%		Concatenate the subsonic and supersonic portions of the solution
		eta = [eta, supsol.x];
		I_1 = [I_1, supsol.y(1,:)];
		I_2 = [I_2, supsol.y(2,:)];
		I_3 = [I_3, supsol.y(3,:)];
		I_4 = [I_4, supsol.y(4,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		SHOCKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	elseif (shocked)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%			First Pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		First pass (inlet -> throat)
%		Do something about BCs		
		if 		(runtype == 1) w_l = [1; NaN; 0; 0]; w_r = [NaN; NaN; NaN; NaN];
		elseif 	(runtype == 2) error('Imposing upstream propagating acoustic wave at the inlet!');
		elseif 	(runtype == 3) w_l = [0; NaN; 1; 0]; w_r = [NaN; NaN; NaN; NaN];	
		elseif	(runtype == 4) w_l = [0; NaN; 0; 1]; w_r = [NaN; NaN; NaN; NaN];
		end

		if (flaggo == 0) 
			error('The splined approach is inappropriate for discontinuous base flows');
		elseif (flaggo == 1)
%			Compute the eta bounds (see appendix of Duran & Moreau) for the supersonic problem
			x_a 	= sqrt(gp1o2*M_a*M_a/(1+gm1o2*M_a*M_a));
			x_bm 	= sqrt(gp1o2*M_b*M_b/(1+gm1o2*M_b*M_b));
			M_p		= sqrt((gm1*M_b*M_b + 2)/(2*gamma*M_b*M_b - gm1));
			x_bp 	= sqrt(gp1o2*M_p*M_p/(1+gm1o2*M_p*M_p));
			x_c		= sqrt(gp1o2*M_c*M_c/(1+gm1o2*M_c*M_c));
			L_1 	= (x_bm - x_a);
			L_2		= (x_bp - x_c);
			L 		= L_1 + L_2;

			x_1		= x_a/L;
			x_2		= L_1/L
			etabounds = [0, x_1-epsilon]
			SPLINES = [M_a M_b M_c]
		else
			error('Need to address other flaggo cases here');
		end

		ODECaller = @(eta, I) DuranMoreauODE(eta, I, Omega, flaggo, SPLINES);
		BCCaller = @(I_a, I_r) ChokedBCs(I_a, I_r, w_l, w_r, etabounds(1), etabounds(2), flaggo, SPLINES);

%		Initialize the BVP
		if (~exist('subsol', 'var'))
			subsol 	= bvpinit(linspace(etabounds(1), etabounds(2), N2), [0.5 0.5 0.5 0.5]);
		end
%		Solve the BVP
		subsol 	= bvp4c(ODECaller, BCCaller, subsol, options);

%		Unpack the BVP solution
		eta = subsol.x;			%spatial coordinate
		I_1 = subsol.y(1,:);	%mass
		I_2 = subsol.y(2,:);	%enthalpy
		I_3 = subsol.y(3,:);	%entropy
		I_4 = subsol.y(4,:);	%composition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%			Second Pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		Do something about the BCs
		I_end 	= subsol.y(:,end);
		[w_l] = charsFromInvrnts(eta(end), I_end, flaggo, SPLINES);
		w_r		= [NaN; NaN; NaN; NaN];

%		Compute the eta bounds (see appendix of Duran & Moreau) for the supersonic problem
		if (flaggo == 0) 
			error('The splined approach is inappropriate for discontinuous base flows');
		elseif  (flaggo == 1)
			etabounds = [x_1+epsilon, x_2-epsilon]
		else
			error('Need to address other flaggo cases here');
		end

		ODECaller = @(eta, I) DuranMoreauODE(eta, I, Omega, flaggo, SPLINES);
		BCCaller = @(I_a, I_r) SupersonicBCs(I_a, I_r, w_l, w_r, etabounds(1), etabounds(2), flaggo, SPLINES);

%		Initialize
		if (~exist('supsol', 'var'))
			supsol 	= bvpinit(linspace(etabounds(1), etabounds(2), N2), [mean(I_1) mean(I_2) mean(I_3) mean(I_4)]);
		end
		supsol 	= bvp4c(ODECaller, BCCaller, supsol, options);

%		Concatenate the subsonic and supersonic portions of the solution
		eta = [eta, supsol.x];
		I_1 = [I_1, supsol.y(1,:)];
		I_2 = [I_2, supsol.y(2,:)];
		I_3 = [I_3, supsol.y(3,:)];
		I_4 = [I_4, supsol.y(4,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%			Third Pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		Do something about the BCs
		I_end 	= subsol.y(:,end);
		w_l 	= charsFromInvrnts(eta(end), I_end, flaggo, SPLINES)
		w_r		= [NaN; NaN; NaN; NaN];

%		Compute the eta bounds (see appendix of Duran & Moreau) for the supersonic problem
		if (flaggo == 0) 
			error('The splined approach is inappropriate for discontinuous base flows');
		elseif(flaggo == 1)
			etabounds = [x_2+epsilon, 1]
		else
			error('Need to address other flaggo cases here');
		end

		ODECaller = @(eta, I) DuranMoreauODE(eta, I, Omega, flaggo, SPLINES);
		BCCaller = @(I_a, I_r) ShockedBCs(I_a, I_r, w_l, w_r, etabounds(1), etabounds(2), flaggo, SPLINES);

%		Initialize
		if (~exist('shksol', 'var'))
			shksol 	= bvpinit(linspace(etabounds(1), etabounds(2), N2), [mean(I_1) mean(I_2) mean(I_3) mean(I_4)]);
		end
		shksol 	= bvp4c(ODECaller, BCCaller, shksol, options);

%		Concatenate the subsonic and supersonic portions of the solution
		eta = [eta, shksol.x];
		I_1 = [I_1, shksol.y(1,:)];
		I_2 = [I_2, shksol.y(2,:)];
		I_3 = [I_3, shksol.y(3,:)];
		I_4 = [I_4, shksol.y(4,:)];
	end

%	Now that the full solution has been obtained, unpack it for plotting
	for i = 1:length(eta)
		I = [I_1(i); I_2(i); I_3(i); I_4(i)];
		[w,s] = charsFromInvrnts(eta(i), I, flaggo, SPLINES);
		if (flaggo == 0)
			EM(i) 	= ppval(SPLINES(1), eta(i));
			psi(i)	= ppval(SPLINES(4), eta(i));
		elseif (flaggo == 1)
			EM(i)	= MFromEtaLVG(eta(i), SPLINES);
			psi(i)	= BaseFlowFromMLVG(EM(i));
		else
			error('You should really address this MfromA stuff');
		end

		phi(i) 		= s(1);
		nu(i) 		= s(2);
		sigma(i)	= s(3);
		xi(i) 		= s(4);

		w_p(i) = w(1);
		w_m(i) = w(2);
		w_s(i) = w(3);
		w_z(i) = w(4);
	end
	disp('DuranMoreau');
%	[x_a x_bm x_bp x_c 0 x_2 1 L_1 L_2 (L_1 + L_2)]
%	[M_a M_b M_p M_c gamma]


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

	if (plot_psi)
		figure();
		plot(eta, psi, 'b', 'linewidth', 2);
		xlabel('$\eta$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		ylabel('$\Psi$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
		set(gca, 'FontSize', 14, 'FontName', 'Times');
	end%(plot_psi)

%	Compute the transfer matrix by storing the values at each end of the solution into the matrix transfer
	transfer = [w_p(1) w_p(end); w_m(1) w_m(end); w_s(1) w_s(end); w_z(1) w_z(end)];
end

