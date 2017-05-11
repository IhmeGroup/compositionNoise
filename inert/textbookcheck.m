function[] = textbookcheck(a, A, MW, Hover)
	close all;

	hh = figure();
	set(hh, 'Position', [0 0 2000 1000]);
	p = 1.0E6;%1 bar in dynes/cm^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	O2 TEST CASE
	T = [0 100 200 298 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800 4000 4400 4800 5200 5600 6000]';
	deltah 	= [-8683 -5777 -2868 0 54 3027 6086 9245 12499 15836 19241 22703 26212 29761 33345 36958 40600 44267 47959 51674 55414 59176 66770 74453 82225 90080 98013 106022 114101 122245 130447 138705 155374 172240 189312 206618 224210]';
	deltah0 = 0;
	MW_BS = 31.999;
	s0		= [0 173.308 193.483 205.148 205.329 213.873 220.693 226.450 231.465 235.920 239.931 243.579 246.923 250.011 252.878 255.556 258.068 260.434 262.673 264.797 266.819 268.748 272.366 275.708 278.818 281.729 284.466 287.050 289.499 291.826 294.043 296.161 300.133 303.801 307.217 310.423 313.457]';
	h_BS = (deltah + deltah0)/MW_BS;
	s0_BS = s0/MW_BS;

%	Generate my solver data (_me)
	h_me = zeros(length(T),1);
	s0_me = zeros(length(T),1);
	Y = [0 0 0 0 0 0 0 1 0];	
	for i = 1:length(T)
	%	[cp, h_me(i), s0_me(i)] = returnSpeciesProperties(T(i), p, Y, a, A, MW, Hover);
		[cp_1(i), h_me(i), s0_me(i), g] = returnSpeciesProperties(T(i), p, Y, a, A, MW)
	end
	psi_1 = (h_me - T.*s0_me)./(cp_1'.*T);
	psi_1 = (h_me)./(cp_1'.*T);

%	Plot the comparison
	subplot(2,3,1);
	plot(T, h_me.*1E-7, 'bp', 'LineWidth', 2);
	hold on;
	plot(T, h_BS, 'k-','LineWidth',2);
	title('$O_2$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$T$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$h$ [kJ/kg]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName','Times');
	aa = legend('My Solver', 'Borgnakke & Sonntag');
	set(aa, 'FontSize', 14, 'FontName', 'Times', 'Location', 'SouthEast');

	subplot(2,3,4);
	plot(T, s0_me.*1E-7, 'bp', 'LineWidth', 2);
	hold on;
	plot(T, s0_BS, 'k-','LineWidth',2);
	title('$O_2$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$T$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$s^o$ [kJ/kg K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName','Times');
	aa = legend('My Solver', 'Borgnakke & Sonntag');
	set(aa, 'FontSize', 14, 'FontName', 'Times', 'Location', 'SouthEast');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	CO2 TEST CASE

%	Generate Borgnakke and Sonntag data (_BS)
	deltah = [-9364 -6457 -3413 0 69 4003 8305 12906 17754 22806 28030 33397 38885 44473 50148 55895 61705 67569 73480 79432 85420 91439 103562 115779 128074 140435 152853 165321 177836 190394 202990 215624 240992 266488 292112 317870 343782]';
	deltah0 = -393522;
	MW_BS = 44.01;
	s0	= [0 179.010 199.976 213.794 214.024 225.314 234.902 243.284 250.752 257.496 263.646 269.299 274.528 279.390 283.931 288.190 292.199 295.984 299.567 302.969 306.207 309.294 315.070 320.384 325.307 329.887 334.170 338.194 341.988 345.576 348.981 352.221 358.266 363.812 368.939 373.711 378.180]';
	h_BS = (deltah + deltah0)/MW_BS;
	s0_BS = s0/MW_BS;

%	Generate my solver data (_me)
	h_me = zeros(length(T),1);
	s0_me = zeros(length(T),1);
	Y = [0 0 1 0 0 0 0 0 0];	
	for i = 1:length(T)
%		[cp, h_me(i), s0_me(i)] = returnSpeciesProperties(T(i), p, Y, a, A, MW, Hover);
		[cp_2(i), h_me(i), s0_me(i), g] = returnSpeciesProperties(T(i), p, Y, a, A, MW)
	end
	psi_2 = (h_me - T.*s0_me)./(cp_2'.*T);
	psi_2 = (h_me)./(cp_2'.*T);

%	Plot the comparison
	subplot(2,3,2);
	plot(T, h_me.*1E-7, 'bp', 'LineWidth', 2);
	hold on;
	plot(T, h_BS, 'k-','LineWidth',2);
	title('$CO_2$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$T$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$h$ [kJ/kg]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName','Times');
	aa = legend('My Solver', 'Borgnakke & Sonntag');
	set(aa, 'FontSize', 14, 'FontName', 'Times', 'Location', 'SouthEast');

	subplot(2,3,5);
	plot(T, s0_me.*1E-7, 'bp', 'LineWidth', 2);
	hold on;
	plot(T, s0_BS, 'k-','LineWidth',2);
	title('$CO_2$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$T$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$s^o$ [kJ/kg K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName','Times');
	aa = legend('My Solver', 'Borgnakke & Sonntag');
	set(aa, 'FontSize', 14, 'FontName', 'Times', 'Location', 'SouthEast');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	H2O TEST CASE

%	Generate Borgnakke and Sonntag data (_BS)
	deltah = [-9904 -6617 -3282 0 62 3450 6922 10499 14190 18002 21937 26000 30190 34506 38941 43491 48149 52907 57757 62693 67706 72788 83153 93741 104520 115463 126548 137756 149073 160484 171981 183552 206892 230456 254216 278161 302295];
	deltah0 = -241826;
	MW_BS = 18.015;
	s0	= [0 152.386 175.488 188.835 189.043 198.787 206.532 213.051 218.739 223.826 228.460 232.739 236.732 240.485 244.035 247.406 250.620 253.690 256.631 259.452 262.162 264.769 269.706 274.312 278.625 282.680 286.504 290.120 293.550 296.812 299.919 302.887 308.448 313.573 318.328 322.764 326.926];

	h_BS = (deltah + deltah0)/MW_BS;
	s0_BS = s0/MW_BS;

%	Generate my solver data (_me)
	h_me = zeros(length(T),1);
	s0_me = zeros(length(T),1);
	Y = [0 0 0 0 0 0 0 0 1];	
	for i = 1:length(T)
%		[cp_3(i), h_me(i), s0_me(i)] = returnSpeciesProperties(T(i), p, Y, a, A, MW, Hover);
		[cp_3(i), h_me(i), s0_me(i), g] = returnSpeciesProperties(T(i), p, Y, a, A, MW)
	end

%	Plot the comparison
	subplot(2,3,3);
	plot(T, h_me.*1E-7, 'bp','LineWidth',2);
	hold on;
	plot(T, h_BS, 'k-','LineWidth',2);
	title('$H_2O$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$T$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$h$ [kJ/kg]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName','Times');
	aa = legend('My Solver', 'Borgnakke & Sonntag');
	set(aa, 'FontSize', 14, 'FontName', 'Times', 'Location', 'SouthEast');

	subplot(2,3,6);
	plot(T, s0_me.*1E-7, 'bp','LineWidth',2);
	hold on;
	plot(T, s0_BS, 'k-','LineWidth',2);
	title('$H_2O$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$T$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$s^o$ [kJ/kg K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName','Times');
	aa = legend('My Solver', 'Borgnakke & Sonntag');
	set(aa, 'FontSize', 14, 'FontName', 'Times', 'Location', 'SouthEast');

	figure();
	plot(T, [cp_1; cp_2; cp_3].*1E-7, 'LineWidth', 3);
	legend('O2', 'CO2', 'H2O');

	figure();
	psi_3 = (h_me-T.*s0_me)./(cp_3'.*T);
	psi_3 = (h_me)./(cp_3'.*T);
	plot(T, [psi_1, psi_2, psi_3], 'LineWidth', 3);
	legend('O2', 'CO2', 'H2O');

end
