function[psi, gctr] = returnPsi(Tbar, pbar, Zbar)
%	Returns the normalized chemical potential and (gibbs free energy, if desired) as a function of mean temperature, pressure, and composition
	global data fuel;
%	Load species properties
	if (fuel == 1)
		data = load('../data/lowStrain/lowStrain.C12H26');
		[Nspecies, species, a, A, MW] = speciesPropsC12H26();
	elseif (fuel == 2)
		data = load('../data/lowStrain/lowStrain.CH4');
		[Nspecies, species, a, A, MW] = speciesPropsCH4();
	elseif (fuel == 3)
		data = load('../data/lowStrain/lowStrain.H2');
		[Nspecies, species, a, A, MW] = speciesPropsH2();
	else
		error(strcat('Undefined fuel type = ', num2str(fuel)));
	end

%	compute the size of the flamelet stored in data
	[N,M] = size(data);

%	Find the Z value to the right of Zbar
	index = 1;
	while (data(1,index+1) < Zbar) %&& ((index + 1) < M) 
		index = index + 1;
		if (index > M)
			break;
		end
	end

%	Obtain the composition on either side of Zbar
	Ylft = data(3:N,index);
	Yrgt = data(3:N,index+1);
	
%	Obtain the mixtur fractions on either side of Zbar
	Zlft = data(1,index);
	Zrgt = data(1,index+1);

%	Central difference to obtain the composition gradient at Zbar
	Yctr = (Zrgt - Zbar)/(Zrgt - Zlft)*Ylft + (Zbar - Zlft)/(Zrgt - Zlft)*Yrgt;

%	Obtain the gibbs free energies at either side of Zbar, and the specific heat at Zbar
	[~, ~, ~, glft] = returnSpeciesProperties(Tbar, pbar, Ylft, a, A, MW);
	[cp, ~, ~, gctr]	= returnSpeciesProperties(Tbar, pbar, Yctr, a, A, MW);
	[~, ~, ~, grgt] = returnSpeciesProperties(Tbar, pbar, Yrgt, a, A, MW);

%	Compute the gradient of the Gibbs free energy at Zbar
	dGdZ = (grgt - glft)/(Zrgt - Zlft);

%	Return the value of psi
	psi = 1/(cp*Tbar)*dGdZ;
end
