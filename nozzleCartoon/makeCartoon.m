function[] = makeCartoon()
	close all;
	h = figure();
	set(h, 'position', [0 0 600, 300]);
%	Generate the nozzle profile
	M_a = 0.3;
	M_b = 2.0;
	N = 101;
	delta = 0.2;
	eps = 0.1;
	gamma = 1.4;
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;
	xihat = linspace(gp1o2*M_a^2/(1+gm1o2*M_a^2), gp1o2*M_b^2/(1+gm1o2*M_b^2), N);
	M = sqrt(2/gp1*xihat./(1-gm1/gp1.*xihat));
	A = (gp1o2).^(-gp1o2/gm1)*(1+gm1o2.*M.*M).^(gp1o2/gm1)./M;
	xi = sqrt(xihat);
	xi = xi - min(xi');
	xi = xi./max(xi);

%	Plot the nozzle
	plot(xi, A/2, 'k-', 'LineWidth', 3);
	hold on;
	plot(xi, -A/2, 'k-', 'LineWidth', 3);	
%	And some upstream piping
	plot([-1.25,0], [1, 1], 'k-', 'LineWidth', 3);
	plot([-1.25,0], -[1, 1], 'k-', 'LineWidth', 3);

%	Add a shock in the downstream portion of the nozzle
	Ny = 51;
	x = xi(64)*ones(Ny,1);
	x(2:Ny-1) = x(2:Ny-1) + eps*rand(Ny-2,1);
	y = linspace(A(64)/2, -A(64)/2, Ny)';
	plot(x, y, 'k--', 'LineWidth', 3);
	
%	Plot some waves entering
	plot([-.5-0*delta, -.5-0*delta], [-1, 1], 'k--', 'LineWidth', 1);
	plot([-.5-1*delta, -.5-1*delta], [-1, 1], 'k--', 'LineWidth', 1);
	plot([-.5-2*delta, -.5-2*delta], [-1, 1], 'k--', 'LineWidth', 1);

%	Plot some waves coming out
	plot([1.5+0*delta, 1.5+0*delta], A(N)/2*[-1, 1], 'k--', 'LineWidth', 1);
	plot([1.5+1*delta, 1.5+1*delta], A(N)/2*[-1, 1], 'k--', 'LineWidth', 1);
	plot([1.5+2*delta, 1.5+2*delta], A(N)/2*[-1, 1], 'k--', 'LineWidth', 1);

%	Plot arrows indicating wave direction
	plot([-.4, -.2], [0, 0], 'k-');
	plot([-.2, -.25], [0, 0.05], 'k-');
	plot([-.2, -.25], [0, -0.05], 'k-');

%	
	plot([2, 2.2], [0, 0], 'k-');
	plot([2.2, 2.15], [0, 0.05], 'k-');
	plot([2.2, 2.15], [0, -0.05], 'k-');

%	Put a box around the flamelet
	plot([-2.75, -1.25], [-0.75, -0.75], 'k-');
	plot([-2.75, -1.25], [0.75, 0.75], 'k-');
	plot([-1.25, -1.25], [0.75, -0.75], 'k-');
	plot([-2.75, -2.75], [0.75, -0.75], 'k-');

%	Plot the T-Z curve to represent the flamelet
	data = load('../lowStrain/lowStrain.CH4')';
	xnorm = data(:,1);
	ynorm = data(:,2);
	xnorm = xnorm - min(xnorm);
	ynorm = ynorm - min(ynorm);
	xnorm = xnorm./max(xnorm);
	ynorm = ynorm./max(ynorm);
	xnorm = xnorm*1.5 - 2.75;
	ynorm  = ynorm*1.5 - 0.75;
	plot(xnorm, ynorm);
	plot(xnorm(150), ynorm(150), 'o');

%	Formatting crap
	set(0, 'defaultTextInterpreter', 'LaTeX');

%	Label the entry and exit states
	text(0, -1.25, '$a$', 'FontSize', 14)
	text(xi(64), -1.25, '$b$', 'FontSize', 14)
	text(1, -1.25, '$c$', 'FontSize', 14)

%	Label the entry and exit perturbation types
	text(-.5-delta, -1.25, '$\xi_a$', 'FontSize', 14);
	text(1.5+delta, -1.25, '$\pi_c^+$', 'FontSize', 14);

%	Label the flamelet
	text(-2.25, 0.5, '$\bar{p}_0, \bar{T}_0$','FontSize', 14);
	text(-2.5, 0.9, 'Flamelet','FontSize', 14);
	text(-2.1, -0.9, '$Z$', 'FontSize', 14);	
	text(-3, 0, '$T$', 'FontSize', 14);	
	text(-1.8, -0.2, '$\bar{Z}$', 'FontSize', 14);

%	Formatting
	axis equal;
	axis([-3,3, -1.5, 1.5])
	axis off;
end
