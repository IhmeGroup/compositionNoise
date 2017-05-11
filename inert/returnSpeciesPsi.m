function[psii] = returnSpeciesPsi(T, p, Ylft, Yctr, Yrgt, Zlft, Zctr, Zrgt, a, A, MW)
	Nspecies = length(MW);
	psii = zeros(Nspecies,1);
	for i = 1:Nspecies
		Y = zeros(Nspecies, 1);
		Y(i) = 1.0;
		[~, ~, ~, g] = returnSpeciesProperties(T, p, Y, a, A, MW);
		g = g.*1E-7;%There are 1E-7 kJ/kg in 1 erg/gram
		dYidZ = (Yrgt(i) - Ylft(i))/(Zrgt - Zlft);
		psii(i) = g*dYidZ;
%		psii(i) = dYidZ;
	end
end
