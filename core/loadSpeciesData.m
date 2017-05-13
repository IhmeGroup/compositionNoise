function[] = loadSpeciesData()
	global mechanism;
	global data;
	global Nspecies;
	global species;
	global a;
	global A;
	global MW;
	global Hover;
	if (mechanism == 1)
		data = load('../data/lowStrain/lowStrain.C12H26');
		[Nspecies, species, a, A, MW, Hover] = speciesPropsC12H26();
	elseif (mechanism == 2)
		data = load('../data/lowStrain/lowStrain.CH4');
		[Nspecies, species, a, A, MW, Hover] = speciesPropsCH4();
	elseif (mechanism == 3)
		data = load('../data/lowStrain/lowStrain.H2');
		[Nspecies, species, a, A, MW, Hover] = speciesPropsH2();
	elseif (mechanism == 4)
		[Nspecies, species, a, A, MW, Hover] = speciesPropsInert();
	elseif (mechanism == 5)%ICSV
		data = load('../data/lowStrain/lowStrain.CH4');
		[Nspecies, species, a, A, MW, Hover] = speciesPropsMW();
	else
		error(strcat('Undefined mechanism type = ', num2str(mechanism)));
	end
end
