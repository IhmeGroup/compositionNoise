function[] = inertDriver()
	close all;
	textbook_check = false;
	sweep_inert = false;
	sweep_H2 = false;
	sweep_CH4 = false;
	sweep_C12H26 = false;
	plot_inert_transfer = false;
	plot_temperature_sensitivity = false;
	plot_pressure_sensitivity = false;
	plot_Scurves = false;
	plot_flamelets = true;
	plot_fuels = false;
	plot_psi = false;
	
%	If desired, plot the validation case of CO2 and H2O enthalpy and entropy vs. Temperature for this program and Borgnakke and Sonntag
	if (textbook_check)
		[Nsp, species, a, A, MW] = speciesPropsInert();
		textbookcheck(a, A, MW);
	end

	if (plot_inert_transfer)
		plotInertTransferFunction();
	end%(plot_inert_transfer)

	if (plot_temperature_sensitivity)
		plotTemperatureSensitivity();
	end%(plot_temperature_sensitivity)

	if (plot_pressure_sensitivity)
		plotPressureSensitivity();
	end%(plot_pressure_sensitivity)

%	If desired, plot the enthalpy and entropy profiles of each species individually. This should be enough to recognize typos in the inputs, since the low temp values should match the high temp values
%	Inert sweep
	if (sweep_inert)
		[Nsp, species, a, A, MW] = speciesPropsInertH2();
		for i = 1:Nsp
			testOneSpecies(i, a, A, MW, species(i))
		end
	end	
%	H2 sweep
	if (sweep_H2);
		[Nsp, species, a, A, MW] = speciesPropsH2();
		for i = 1:Nsp
			testOneSpecies(i, a, A, MW, species(i))
		end
	end%(sweep_H2)
%	CH4 sweep
	if (sweep_CH4)
		[Nsp, species, a, A, MW] = speciesPropsCH4();
		for i = 1:10
			testOneSpecies(i, a, A, MW, species(i));
		end
	end%(sweep_CH4)
%	C12H26 sweep
	if (sweep_C12H26)
		[Nsp, species, a, A, MW] = speciesPropsC12H26();
		for i = 1:10
			testOneSpecies(i, a, A, MW, species(i));
		end
	end

	if (plot_Scurves)
		plotSCurves();
	end%(plot_Scurves)

	if (plot_flamelets)
		data = load('./lowStrain/lowStrain.H2');
		plotH2Flamelet(data);
%		print -depsc H2Flamelet.eps
		data = load('./lowStrain/lowStrain.CH4');
		plotCH4Flamelet(data);
%		print -depsc CH4Flamelet.eps
		data = load('./lowStrain/lowStrain.C12H26');
		plotC12H26Flamelet(data);
%		print -depsc C12H26Flamelet.eps
	end%(plot_flamelets)

	if (plot_fuels)
		plotFuelTransferFunctions();
	end

	if (plot_psi)
		plotFuelPsiFields();
	end

end%inertDriver()
