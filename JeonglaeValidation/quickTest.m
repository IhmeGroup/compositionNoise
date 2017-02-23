function[] = quickTest()
	close all;
	global fuel;
	fuel = 4;
	[gamma, T0, p0, Zbar] = returnAmbientState();

	for i = 100:-1:1
		T(i) = i/100*T0;
		psi(i) = returnPsi(T(i), 101000, 0.5);
	end

	plot(T, psi, 'LineWidth', 4);
	xlabel('T [K]');
	ylabel('$\Psi$', 'Interpreter', 'LaTeX');
	set(gca, 'FontSize', 14);
end
