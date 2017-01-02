function[data] = loadFuelData(fuel)
	if (fuel == 1) % Dodecane
		data = load('../data/lowStrain/lowStrain.C12H26');
	elseif (fuel == 2) % Methane
		data = load('../data/lowStrain/lowStrain.CH4');
	elseif (fuel == 3) % 
		data = load('../data/lowStrain/lowStrain.H2');
	end
end
