function[] = generateHelmholtzSweepData(fueltype);
%	Plots the spatial evolution of the Duran and Moreau solution for linear velocity mach profiles at a variety of Helmholtz numbers
%	The degree of oscillation can be (poorly) controlled by changing the value for epsilon. For the case considered here, 1E-8 seems about optimal. 
%	If you change this parameter, then you're also going to have to change the tolerances on the wrapping function and possibly the ylim as well.
%	Sweeping over the dimensionless frequency in increments of 0.5

%	If fuel == 1, C12H26
%	If Fuel == 2, CH4
%	If Fuel == 3, H2

	global fuel;
	fuel = fueltype;

	addpath('../core');
	addpath('../data');

	WP = cell(5,1);
	WM = cell(5,1);
	WS = cell(5,1);
	WZ = cell(5,1);
	ETA = cell(5,1);


	for i = 1:5
		Omega = (i-1)*0.5;
%		Solve the ODE and return the solution in terms of the characteristic variables
		if (i == 1)
			[~, subsol, supsol, eta, w_p, w_m, w_s, w_z, SPLINES] = DuranMoreau(0.29, 1.5, 0, Omega, 4, true);
		else
			[~, subsol, supsol, eta, w_p, w_m, w_s, w_z, SPLINES] = DuranMoreau(0.29, 1.5, 0, Omega, 4, true, SPLINES, subsol, supsol);
		end
		WP(i) 	= {w_p};
		WM(i) 	= {w_m};
		WS(i) 	= {w_s};
		WZ(i) 	= {w_z};
		ETA(i) 	= {eta};
	end

	save('HelmholtzSweepData.mat', 'ETA', 'WP', 'WM', 'WS', 'WZ');

end
