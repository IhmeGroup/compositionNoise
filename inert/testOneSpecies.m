function[] = testOneSpecies(index, a, A, MW, name)
	p = 1.0E6;%1 bar in dynes/cm^2
	T = [0 100 200 298 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800 4000 4400 4800 5200 5600 6000]';

%	Generate my solver data (_me)
	h_me = zeros(length(T),1);
	s0_me = zeros(length(T),1);
	Y = [0 0 0 0 0 0 0 0 0];	
	Y(index) = 1;
	for i = 1:length(T)
		[cp, h_me(i), s0_me(i)] = returnSpeciesProperties(T(i), p, Y, a, A, MW); 
	end

%	Plot the comparison
	hh = figure();
	set(hh, 'Position', [0 0 1000 400]);
	subplot(1,2,1);
	plot(T, h_me.*1E-7, 'bp', 'LineWidth', 2);
	title(strcat('$', name, '$'), 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$T$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$h$ [kJ/kg]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName','Times');

	subplot(1,2,2);
	plot(T, s0_me.*1E-7, 'bp', 'LineWidth', 2);
	title(strcat('$',name,'$'), 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$T$ [K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	ylabel('$s^o$ [kJ/kg K]', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	set(gca, 'FontSize', 14, 'FontName','Times');
end
