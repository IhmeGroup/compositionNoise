function[] = plotFlamelet()
	global fuel;
	close all;
	fs = 22;
	lw = 4;
%	Plot the flamelet structure figure
	for fuel = 1:3
		if (fuel == 1)%C12H26
			data2 = load('../data/lowStrain/lowStrain.C12H26');
		elseif (fuel == 2)%CH$
			data2 = load('../data/lowStrain/lowStrain.CH4');
		elseif (fuel == 3)%H2
			data2 = load('../data/lowStrain/lowStrain.H2');
		end
		data2 = data2';
		Z = data2(:,1);
		T = data2(:,2);

		if (fuel == 1)%C12H26
			Y_fuel = data2(:,24);
			Y_H2O = data2(:,8);
			Y_O2 = data2(:,10);
			Z_st = 0.0627964;
		elseif (fuel == 2)%CH4
			Y_fuel = data2(:,18);
			Y_H2O = data2(:,31);
			Y_O2 = data2(:,5);
			Z_st = 0.0551538;
		elseif (fuel == 3)%H2
			Y_fuel = data2(:,9);
			Y_H2O = data2(:,10);
			Y_O2 = data2(:,6);
			Z_st = 0.0285207;
		end
		phi = Z./(1 - Z).*(1-Z_st)./Z_st;
		x = phi./(1+phi);
		hh = figure();
		set(hh, 'Position', [0 0 650 450]);
		[ax, h1, h2] = plotyy(x, T, x, [Y_fuel, Y_O2, Y_H2O]);
		hold on;
		plot([0.5 0.5], [0, 2500], 'k--','LineWidth', lw);
		set(h1, 'LineWidth', lw, 'LineStyle', '-','Color','k');
		set(h2(1), 'LineWidth', lw, 'LineStyle', '--','color', 'b');
		set(h2(2), 'LineWidth', lw, 'LineStyle', '-.','color','r');
		set(h2(3), 'LineWidth', lw, 'LineStyle', ':','color','m');
		xlabel('$Z^* [-]$', 'FontSize', fs, 'FontName', 'Times', 'Interpreter','LaTeX');
		ylabel(ax(1), '$T$ [K]','FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
		ylabel(ax(2), '$Y_i [-]$','FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
		set(ax(1), 'FontSize', fs, 'FontName', 'Times','YColor', 'k');
		set(ax(2), 'FontSize', fs, 'FontName', 'Times', 'YColor', 'k');
		if (fuel == 1)
			print -depsc C12H26FlameStructure.eps
		elseif (fuel == 2)
			print -depsc CH4FlameStructure.eps
		elseif (fuel == 3)
			print -depsc H2FlameStructure.eps
		end
	end%for fuel = 1:3
end
