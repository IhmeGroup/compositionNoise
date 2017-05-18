function[OMEGA TRANS PHASE] = generateJeonglaeSweepData(forcing)
%	This function is intended to validate the (non-composition) portion of the solver by replicating the data found in
%	Duran and Moreau, Solution of the quasi-one-dimensional linearized Euler equations using flow invariants and the Magnus expansion, JFM vol 723, pp 190-231. 2013
%	The figure number (i.e. 5a, 5b, etc.) correspond to the equivalent figure in Duran and Moreau
%	After generating the equivalent data using O'Brien's solver, results pulled from D&M's paper are then superimposed on the plot at the end of the code
%	This data was taken using grabit in MATLAB and is thus fairly imprecise. Discrepencies between D&M and this solver can be attributed partially to this and partially to numerical error.
%	Varying the ``tightness'' of the throat perturbation parameter epsilon is found to have an appreciable effect on the agreement (small values of epsilon are better), but this also slows down the code
%	Also, a subsonic case is included b/c it seemed like a valuable addition.

	addpath('../core');
	addpath('../data');
	addpath('../speciesProps');

%	parameters
	suppress = true;
	mechanism = 6;
	flaggo = 1;
	Nomega = 51;
	dOmega = 2/(Nomega-1);

%	Allocate data for efficiency
	TRANS = zeros(3,Nomega,5);
	PHASE = zeros(3, Nomega, 3);

	for test = 1
		M_a = 1.1;
		M_b = 1.5;
		M_c = 0.0;

%	for test = 1:3
%		if (test == 1)		M_a = 0.29; M_b = 0.88; M_c = 0.0;
%		elseif (test == 2) 	M_a = 0.29; M_b = 1.02; M_c = 0.0; 
%		elseif (test == 3) 	M_a = 0.29; M_b = 1.50; M_c = 0.0; 
%		end

%		Build the base flow
		count = 0;
		for omega = 0:dOmega:2.5
			disp(omega)
			count = count + 1;
			OMEGA(count) = omega;
			if (exist('subsol') && exist('supsol') && exist('shksol'))
				[transfer, subsol, supsol, shksol, SPLINES] = DuranMoreau(M_a, M_b, M_c, omega, forcing, mechanism, suppress, flaggo, SPLINES, subsol, supsol, shksol);
			else
%					1		2		3		4		5				       1    2    3     4       5         6          7        8	
				[transfer, subsol, supsol, shksol, SPLINES] = DuranMoreau(M_a, M_b, M_c, omega, forcing, mechanism, suppress, flaggo);
				
			end
			TRANS(test,count,:) 	= [transfer(1,2), transfer(2,2), transfer(3,2), transfer(4,2), transfer(2,1)];
			PHASE(test, count, :) 	= [	atan2(imag(transfer(1,2)), real(transfer(1,2))), ...
									 	atan2(imag(transfer(2,2)), real(transfer(2,2))), ...
										atan2(imag(transfer(2,1)), real(transfer(2,1)))];
		end%omega
	end%for test
	save('JeonglaeData.mat', 'OMEGA', 'TRANS', 'PHASE');
end
