function[] = testSpeciesProperties(mechanism)
	N = 1001;
	addpath('../speciesProps/');

	if (mechanism == 1)
		[Nspecies, species, a, A, MW, Hover] = speciesPropsC12H26();
	elseif (mechanism == 2)
		[Nspecies, species, a, A, MW, Hover] = speciesPropsCH4();
	elseif (mechanism == 3)
		[Nspecies, species, a, A, MW, Hover] = speciesPropsH2();
	elseif (mechanism == 4)
		[Nspecies, species, a, A, MW, Hover] = speciesPropsInert();
	elseif (mechanism == 5)
		[Nspecies, species, a, A, MW, Hover] = speciesPropsMW();
	elseif (mechanism == 6)
		[Nspecies, species, a, A, MW, Hover] = speciesPropsNobleGas();
	end

	T = linspace(300,5000,N);
	b = figure();
	for i = 1:Nspecies
		Y = zeros(Nspecies,1);
		Y(i) = 1;
		for j = 1:N
			[cp(j), h(j), s(j), g(j)] = returnSpeciesProperties(T(j), 1E5, Y, a, A, MW, Hover);
		end
		figure(b);
		subplot(2,2,1);
		plot(T, cp)
		xlabel('T [K]');
		ylabel('cp [J/kg*K]');
		subplot(2,2,2);
		plot(T, h);
		xlabel('T [K]');
		ylabel('h [J/kg]');
		subplot(2,2,3)
		plot(T, s);
		xlabel('T [K]');
		ylabel('s [J/kg*K]');
		subplot(2,2,4);
		plot(T, g);
		xlabel('T [K]');
		ylabel('g [J/kg]');

		me = cell2mat((species(i)));
		i
		msg = strcat('Does ', me,  'look OK?')
		input(msg)
	end
	
end
