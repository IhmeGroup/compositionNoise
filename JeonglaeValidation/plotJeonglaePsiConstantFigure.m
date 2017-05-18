function[] = plotJeonglaePsiConstantFigure()
	close all;
	data00 = load('JeonglaeData.PsiConstant.0.mat');
	data01 = load('JeonglaeData.PsiConstant.p1.mat');
	data20 = load('JeonglaeData.PsiConstant.m20.mat');
	data40 = load('JeonglaeData.PsiConstant.m40.mat');
	dataJ	= load('data_LEE.dat');

	figure();
	plot(data00.OMEGA(:), abs(data00.TRANS(1,:,1))./2, 'g-', 'LineWidth', 3);
	hold on;
	plot(data01.OMEGA(:), abs(data01.TRANS(1,:,1))./2, 'b-', 'LineWidth', 3);
	plot(data20.OMEGA(:), abs(data20.TRANS(1,:,1))./2, 'c-', 'LineWidth', 3);
	plot(data40.OMEGA(:), abs(data40.TRANS(1,:,1)./2), 'm-', 'LineWidth', 3);
	hold on;
	dataJ(:,1) = dataJ(:,1);
	plot(dataJ(:,1), abs(dataJ(:,2)), 'gx', 'MarkerSize', 15);
	plot(dataJ(:,1), abs(dataJ(:,3)), 'bx', 'MarkerSize', 15);
	plot(dataJ(:,1), abs(dataJ(:,4)), 'cx', 'MarkerSize', 15);
	plot(dataJ(:,1), abs(dataJ(:,5)), 'mx', 'MarkerSize', 15);

	h = legend('$\Psi = 0$', '$\Psi = 1$', '$\Psi = -20$', '$\Psi = -40$', 'Location', 'NorthEast');
	set(h, 'Interpreter', 'LaTeX', 'FontSize', 16, 'FontName', 'Times');
	xlabel('He', 'FontName', 'Times', 'FontSize', 16);
	ylabel('$\pi^+_b  \slash \xi_a$', 'Interprete', 'LaTeX', 'FontSize', 16, 'FontName', 'Times');
	set(gca, 'FontSize', 16, 'FontName', 'Times');
end
