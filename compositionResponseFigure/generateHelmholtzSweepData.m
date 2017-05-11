function[OMEGA TRANS PHASE] = validateDuranMoreau(forcing)
%	This function is intended to validate the (non-composition) portion of the solver by replicating the data found in
%	Duran and Moreau, Solution of the quasi-one-dimensional linearized Euler equations using flow invariants and the Magnus expansion, JFM vol 723, pp 190-231. 2013
%	The figure number (i.e. 5a, 5b, etc.) correspond to the equivalent figure in Duran and Moreau
%	After generating the equivalent data using O'Brien's solver, results pulled from D&M's paper are then superimposed on the plot at the end of the code
%	This data was taken using grabit in MATLAB and is thus fairly imprecise. Discrepencies between D&M and this solver can be attributed partially to this and partially to numerical error.
%	Varying the ``tightness'' of the throat perturbation parameter epsilon is found to have an appreciable effect on the agreement (small values of epsilon are better), but this also slows down the code
%	Also, a subsonic case is included b/c it seemed like a valuable addition.

	addpath('../core');
	addpath('../data');

%	parameters
	Nomega = 21;
	dOmega = 2/(Nomega-1);
	global fuel data beta;
	beta = -2;%Use linear vel grad geometry
	data = loadFuelData(fuel);

%	Allocate data for efficiency
	TRANS = zeros(3,Nomega,5);
	PHASE = zeros(3, Nomega, 3);

	for test = 1:1
		if (test == 1)		M_a = 1.10; M_b = 1.50; M_c = 0.0;
		elseif (test == 2) 	M_a = 0.29; M_b = 1.02; M_c = 0.0; 
		elseif (test == 3) 	M_a = 0.29; M_b = 1.50; M_c = 0.0; 
		end

		global param;
		param = zeros(17,1);
		[gamma, T0, p0, Zbar] = returnAmbientState();
%				 1	  2    3    4      5  6   7   8     9 10 11 12 13 14 15 16 17 
		param = [M_a; M_b; M_c; gamma; 0; T0; p0; Zbar; 0; 0; 0; 0; 0; 0; 0; 0; 0];
		[SPLINES] = buildBaseFlowSplines()
%		Build the base flow
		count = 0;
		for omega = 0:dOmega:2.0
			disp(omega)
			count = count + 1;
			OMEGA(count) = omega;
%			if ((exist('subsol', 'var')) && (exist('supsol', 'var')))
%				[transfer, subsol, supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, M_c, omega, forcing, true, SPLINES, subsol, supsol);
%			else
				[transfer, subsol, supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, M_c, omega, forcing, true, SPLINES);
%			end
			TRANS(test,count,:) = [transfer(1,2), transfer(2,2), transfer(3,2), transfer(4,2), transfer(2,1)];
			PHASE(test, count, :) = [	atan2(imag(transfer(1,2)), real(transfer(1,2))), ...
									 	atan2(imag(transfer(2,2)), real(transfer(2,2))), ...
										atan2(imag(transfer(2,1)), real(transfer(2,1)))];
		end%omega
		clear subsol supsol;
	end%for test
%	save('validationData.mat', 'OMEGA', 'TRANS', 'PHASE');
end
