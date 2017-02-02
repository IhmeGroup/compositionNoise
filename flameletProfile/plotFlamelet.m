function[] = plotFlamelet(fueltype)
	global fuel;
	fuel = fueltype;
	close all;
%	Plot the flamelet structure figure
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
	elseif (fuel == 2)%CH4
		Y_fuel = data2(:,18);
		Y_H2O = data2(:,31);
		Y_O2 = data2(:,5);
	elseif (fuel == 3)%H2
		Y_fuel = data2(:,9);
		Y_H2O = data2(:,10);
		Y_O2 = data2(:,6);
	end
	Z_st = 0.0627964;
	x = Z./(Z + Z_st);
	hh = figure();
	subplot(1,2,1);
%	set(hh, 'Position', [0 0 650 450]);
	set(hh, 'Position', [0 0 1600 500]);
	[ax, h1, h2] = plotyy(x, T, x, [Y_fuel, Y_O2, Y_H2O]);
	hold on;
	plot([0.5 0.5], [0, 2500], 'k--','LineWidth', 2);
	set(h1, 'LineWidth', 2, 'LineStyle', '-','Color','k');
	set(h2(1), 'LineWidth', 2, 'LineStyle', '--','color', 'b');
	set(h2(2), 'LineWidth', 2, 'LineStyle', '-.','color','r');
	set(h2(3), 'LineWidth', 2, 'LineStyle', ':','color','m');
	xlabel('$Z \slash (Z + Z_{st}) [-]$', 'FontSize', 16, 'FontName', 'Times', 'Interpreter','LaTeX');
	ylabel(ax(1), '$T$ [K]','FontSize', 16, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	ylabel(ax(2), '$Y_i [-]$','FontSize', 16, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	set(ax(1), 'FontSize', 16, 'FontName', 'Times','YColor', 'k');
	set(ax(2), 'FontSize', 16, 'FontName', 'Times', 'YColor', 'k');
%	print -depsc Flamelet.eps

%	Plot the Psi-Gibbs figure, b/c that is now a thing we're including
	addpath('../core');

	N = length(T)
	pbar = 1E5;
	psi = zeros(N,1);
	gctr = zeros(N,1);
	for i = 2:N-1
		[psi(i), gctr(i)] = returnPsi(T(i), pbar, Z(i) + 1E-5);
	end
	psi(1) = (psi(3) - psi(2))/(Z(3) - Z(2))*(-Z(2)) + psi(2);
	psi(N) = (psi(N-1) - psi(N-2))/(Z(N-1) - Z(N-2))*(Z(N)-Z(N-2)) + psi(N-1);
	gctr(1) = (gctr(3) - gctr(2))/(Z(3) - Z(2))*(-Z(2)) + gctr(2);
	gctr(N) = (gctr(N-1) - gctr(N-2))/(Z(N-1) - Z(N-2))*(Z(N)-Z(N-2)) + gctr(N-1);
	
%	gg = figure();
	subplot(1,2,2);
%	set(gg, 'Position', [0 0 650 450]);
	gctr = gctr/1E10;
	[ax, h1, h2] = plotyy(x, psi, x, gctr);
	set(h1, 'LineWidth', 2, 'LineStyle', '-','Color','k');
	set(h2(1), 'LineWidth', 2, 'LineStyle', '--','color', 'b');
%	set(h2(2), 'LineWidth', 2, 'LineStyle', '-.','color','r');
%	set(h2(3), 'LineWidth', 2, 'LineStyle', ':','color','m');
	xlabel('$Z \slash (Z + Z_{st}) [-]$', 'FontSize', 16, 'FontName', 'Times', 'Interpreter','LaTeX');
	ylabel(ax(1), '$\Psi [-]$','FontSize', 16, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	ylabel(ax(2), '$g$ [MJ/kg]','FontSize', 16, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	set(ax(1), 'FontSize', 16, 'FontName', 'Times','YColor', 'k');
	set(ax(2), 'FontSize', 16, 'FontName', 'Times', 'YColor', 'k');
	
	print -depsc Flamelet.eps


end
