function[transfer] = singleRunTest(notbeta, omega)
	addpath('../core');
	addpath('../data');

%	parameters
	global fuel data beta;
	beta = notbeta;
	fuel = 2;
	data = loadFuelData(fuel);

%	Allocate data for efficiency
	M_a = 0.3059;
	M_b = 2.1971;
	M_c = 0;
	forcing = 4;%Dump some composition waves in there

	[gamma, T0, p0, Zbar] = returnAmbientState();

	global param;
	param = zeros(17,1);
%			 1	  2    3    4      5  6   7   8     9 10 11 12 13 14 15 16 17 
	param = [M_a; M_b; M_c; gamma; 0; T0; p0; Zbar; 0; 0; 0; 0; 0; 0; 0; 0; 0];
	[SPLINES] = buildBaseFlowSplines()
%	Build the base flow
	if ((exist('subsol', 'var')) && (exist('supsol', 'var')))
		[transfer, subsol, supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, M_c, omega, forcing, true, SPLINES, subsol, supsol);
	else
		[transfer, subsol, supsol, ~, ~, ~, ~, ~, SPLINES] = DuranMoreau(M_a, M_b, M_c, omega, forcing, true, SPLINES);
	end
	transfer = abs(transfer);
%	TRANS(i,count,:) = [transfer(1,2), transfer(2,2), transfer(3,2), transfer(4,2), transfer(2,1)];
%	PHASE(i, count, :) = [	atan2(imag(transfer(1,2)), real(transfer(1,2))), ...
%						 	atan2(imag(transfer(2,2)), real(transfer(2,2))), ...
%							atan2(imag(transfer(2,1)), real(transfer(2,1)))];
end
