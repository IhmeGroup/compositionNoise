function[] = singleComponentValidation(test)
	close all;
	if (test == 1)
%		Borgnakke & Sonntag			
		T_N2 = [0; 100; 200; 298; 300; 400; 500; 600; 700; 800; 900; 1000; 1100; 1200; 1300; 1400; 1500; 1600; 1700; 1800; 1900; 2000; 2200; 2400; 2600; 2800; 3000; 3200; 3400; 3600; 3800; 4000; 4400; 4800; 5200; 5600; 6000];
		MW_N2 = 28.013;%kg/kmol
		h0_N2 = 0;
		h_N2 = ([-8670; -5768; -2857; 0; 54; 2971; 5911; 8894; 11937; 15046; 18223; 21463; 24760; 28109; 31503; 34936; 38405; 41904;45430; 48979; 52549; 56137; 63362; 70640; 77963; 85323; 92715; 100134; 107577; 115042; 122526; 130027; 145078; 160188; 175352; 190572; 205848] + h0_N2)./MW_N2;
		s_N2 = [0; 159.812; 179.985; 191.609; 191.789; 200.181; 206.740; 212.177; 216.865; 221.016; 224.757; 228.171; 231.314; 234.227; 236.943; 239.487; 241.881; 244.139; 246.276; 248.304; 250.234; 252.075; 255.518; 258.684; 261.615; 264.342; 266.892; 269.286; 271.542; 273.675; 275.698; 277.622; 281.209; 284.495; 287.530; 290.349; 292.984]./MW_N2;
		g_N2 = [h_N2 - T_N2.*s_N2];
		N_N2 = length(T_N2);
		h0_N2 = ones(N_N2,1)*h0_N2;
		[Nspecies, species, a, A, MW, Hover] = speciesPropsC12H26();
		h = zeros(N_N2,1);
		s = zeros(N_N2,1);
		g = zeros(N_N2,1);
%		My data		
		Y_N2 = zeros(Nspecies, 1);
		N2 = 24;
		Y_N2(N2) = 1;
		h0 = ones(N_N2,1)*Hover(N2);
		for i = 1:N_N2
			[~, h(i), s(i), g(i)] = returnSpeciesProperties(T_N2(i), 1.01E5, Y_N2, a, A, MW);
		end
		figure();
		subplot(1,3,1);
		plot(T_N2, h, 'LineWidth', 3);
		hold on;
		plot(T_N2, h0, 'LineWidth', 3);
		plot(T_N2, h_N2, 'rx', 'MarkerSize', 14);
		plot(T_N2, h0_N2, 'rx', 'MarkerSize', 14);
		xlabel('T [K]');
		ylabel('h [kJ/kg]');
		subplot(1,3,2);
		plot(T_N2, s, 'LineWidth', 3);
		hold on;
		plot(T_N2, s_N2, 'rx', 'MarkerSize', 14);
		xlabel('T [K]');
		ylabel('s [kJ/kg*K]');
		subplot(1,3,3);
		plot(T_N2, g, 'LineWidth', 3);
		hold on;
		plot(T_N2, g_N2, 'rx', 'MarkerSize', 14);
		xlabel('T [K]');
		ylabel('g [kJ/kg]');
	elseif (test == 2)
%		Borgnakke & Sonntag			
		T_O2 = [0; 100; 200; 298; 300; 400; 500; 600; 700; 800; 900; 1000; 1100; 1200; 1300; 1400; 1500; 1600; 1700; 1800; 1900; 2000; 2200; 2400; 2600; 2800; 3000; 3200; 3400; 3600; 3800; 4000; 4400; 4800; 5200; 5600; 6000];
		MW_O2 = 31.999;
		h0_O2 = 0;
		h_O2 = ([-8683; -5777; -2868; 0; 54; 3027; 6086; 9245; 12499; 15836; 19241; 22703; 26212; 29761; 33345; 36958; 40600; 44267; 47959; 51674; 55414; 49176; 66770; 74453; 82225; 90080; 98013; 106022; 114101; 122245; 130447; 138705; 155374; 172240; 189312; 206618; 224210]+ h0_O2)./MW_O2;
		s_O2 = [0; 173.308; 193.483; 205.148; 205.329; 213.873; 220.693; 226.450; 231.465; 235.920; 239.931; 243.579; 246.923; 250.011; 252.878; 255.556; 258.068; 260.434; 262.673; 264.797; 266.819; 268.748; 262.366; 275.708; 278.818; 281.729; 284.466; 287.050; 289.499; 291.826; 294.043; 296.161; 300.133; 303.801; 307.217; 310.423; 313.457]./MW_O2;
		g_O2 = [h_O2 - T_O2.*s_O2];
		N_O2= length(T_O2);
		h0_O2 = ones(N_O2,1)*h0_O2;
		[Nspecies, species, a, A, MW, Hover] = speciesPropsC12H26();
		h = zeros(N_O2,1);
		s = zeros(N_O2,1);
		g = zeros(N_O2,1);
%		My data		
		Y_O2 = zeros(Nspecies, 1);
		O2 = 8;
		Y_O2(O2) = 1;
		h0 = ones(N_O2,1)*Hover(O2);
		for i = 1:N_O2
			[~, h(i), s(i), g(i)] = returnSpeciesProperties(T_O2(i), 1.01E5, Y_O2, a, A, MW);
		end
		figure();
		subplot(1,3,1);
		plot(T_O2, h, 'LineWidth', 3);
		hold on;
		plot(T_O2, h0, 'LineWidth', 3);
		plot(T_O2, h_O2, 'rx', 'MarkerSize', 14);
		plot(T_O2, h0_O2, 'rx', 'MarkerSize', 14);
		xlabel('T [K]');
		ylabel('h [kJ/kg]');
		subplot(1,3,2);
		plot(T_O2, s, 'LineWidth', 3);
		hold on;
		plot(T_O2, s_O2, 'rx', 'MarkerSize', 14);
		xlabel('T [K]');
		ylabel('s [kJ/kg*K]');
		subplot(1,3,3);
		plot(T_O2, g, 'LineWidth', 3);
		hold on;
		plot(T_O2, g_O2, 'rx', 'MarkerSize', 14);
		xlabel('T [K]');
		ylabel('g [kJ/kg]');
	elseif (test == 3)
	end
end
