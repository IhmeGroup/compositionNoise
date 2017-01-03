function[] = generateMaMbContourData()
%	Plots the spatial evolution of the Duran and Moreau solution for linear velocity mach profiles at a variety of Helmholtz numbers
%	The degree of oscillation can be (poorly) controlled by changing the value for epsilon. For the case considered here, 1E-8 seems about optimal. 
%	If you change this parameter, then you're also going to have to change the tolerances on the wrapping function and possibly the ylim as well.
	close all;

	addpath('../core');
	addpath('../data');

	Na = 21;
	Nb = 21;
	MAMAX = 1;
	MBMAX = 2;

	WPZ = zeros(Na, Nb);
	WMZ = zeros(Na, Nb);
	WSZ = zeros(Na, Nb);
	WZZ = zeros(Na, Nb);
	WPS = zeros(Na, Nb);
	WMS = zeros(Na, Nb);
	WSS = zeros(Na, Nb);
	WZS = zeros(Na, Nb);
	MA = zeros(Na, Nb);
	MB = zeros(Na, Nb);
	MBEX = zeros(Na, Nb);
	epsilon = 4E-8;

%	parameters
	global fuel data param;
	fuel = 2;
	data = loadFuelData(fuel);

	[gamma, T0, p0, Zbar] = returnAmbientState();

%	Pre-compute for compactness
	gp1 = gamma + 1;
	gm1 = gamma - 1;
	gp1o2 = gp1/2;
	gm1o2 = gm1/2;

	suppress = true;

	B = [1:Nb, Nb:-1:1];

%	Sweeping over the dimensionless frequency 
	for i = 1:5
		Omega = (i-1)*0.5
%		Solve the ODE and return the solution in terms of the characteristic variables
%		Loop over inlet speed
		for a = 1:Na
			M_a = MAMAX/(Na-1)*(a-1)
			if (M_a == 0) M_a = 0.01; end
			if (M_a == 1) M_a = 0.99; end

%			Loop over a snake pattern so that the prevoius solution is always a high quality IC for the next solution
			for bcounter = 1:2*Nb
				b = B(bcounter)
				if (bcounter == Nb+1) a = a + 1; end
				M_b = MBMAX/(Nb-1)*(b-1)
				if (M_b == 1) M_b = 0.99; end
				if (M_b == 0) M_b = 0.01; end

				if (M_b ~= M_a)
					param = zeros(17,1);
%							 1    2    3  4		 5  6   7   8     9 10 11 12 13 14 15 16 17
					param = [M_a; M_b; 0; gamma; 0; T0; p0; Zbar; 0; 0; 0; 0; 0; 0; 0; 0; 0];
					[SPLINES] = buildBaseFlowSplines();
					if (a == 1) && (b == 2)
						[z_transfer, z_subsol, z_supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, 0, Omega, 3, true, SPLINES);
						[s_transfer, s_subsol, s_supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, 0, Omega, 4, true, SPLINES);
					else
						[z_transfer, z_subsol, z_supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, 0, Omega, 3, true, SPLINES, z_subsol, z_supsol);
						[s_transfer, s_subsol, s_supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, 0, Omega, 4, true, SPLINES, s_subsol, s_supsol);
					end
				else
					z_transfer = zeros(4,2);
					z_transfer(1,2) = 0;%w_p = 0
					z_transfer(2,2) = 0;%w_m = 0
					z_transfer(3,2) = 0;%w_s = 0
					z_transfer(4,2) = 0;%w_z = 0
					s_transfer = zeros(4,2);
					s_transfer(1,2) = 0;%w_p = 0
					s_transfer(2,2) = 0;%w_m = 0
					s_transfer(3,2) = 0;%w_s = 0
					s_transfer(4,2) = 0;%w_z = 0
				end
				WPZ(a,b) = z_transfer(1,2);
				WMZ(a,b) = z_transfer(2,2);
				WSZ(a,b) = z_transfer(3,2);
				WZZ(a,b) = z_transfer(4,2);
				WPS(a,b) = s_transfer(1,2);
				WMS(a,b) = s_transfer(2,2);
				WSS(a,b) = s_transfer(3,2);
				WZS(a,b) = s_transfer(4,2);
				MA(a,b)  = M_a;
				MB(a,b)	 = M_b;
			end%for M_b
		end%for M_a

%		The problem is ill-determined for nozzles with M_a == M_b, so linearly interpolate these values from their neighbors
		for i = 2:Na-1
			WPS(i,i) = 0.25*(WPS(i,i+1) + WPS(i,i-1) + WPS(i+1,i) + WPS(i-1,i));
			WPZ(i,i) = 0.25*(WPZ(i,i+1) + WPZ(i,i-1) + WPZ(i+1,i) + WPZ(i-1,i));
		end
		WPS(1,1) = (WPS(1,2) + WPS(2,1))/2;
		WPS(Na, Na) = (WPS(Na,Na-1) + WPS(Na-1, Na))/2;
		WPZ(1,1) = (WPZ(1,2) + WPZ(2,1))/2;
		WPZ(Na, Na) = (WPZ(Na,Na-1) + WPZ(Na-1, Na))/2;

		filename = strcat('MaMbContoursData.', num2str(Omega), '.mat');
		save(filename, 'MA', 'MB', 'WPZ', 'WPS');
	end%for Omega
end
