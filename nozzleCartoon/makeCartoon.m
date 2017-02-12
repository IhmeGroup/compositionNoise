function[] = makeCartoon()
	close all;

	addpath('../data');

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
%	And some downstream piping
	plot([xi(end), xi(end)+1.5], [A(end)/2, A(end)/2], 'k-', 'LineWidth', 3);
	plot([xi(end), xi(end)+1.5], -[A(end)/2, A(end)/2], 'k-', 'LineWidth', 3);

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

%	Label the entry and exit states
	text(0, -1.25, '$a$', 'FontSize', 14)
	text(xi(64), -1.25, '$b$', 'FontSize', 14)
	text(1, -1.25, '$c$', 'FontSize', 14)

%	Label the entry and exit perturbation types
	text(-.5-delta, -1.25, '$\xi_a$', 'FontSize', 14);
	text(1.5+delta, -1.25, '$\pi_c^+$', 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	PLOT THE FLAMELET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Specify some parameters for the flamelet plot
	llcx = -3.00;
	llcy = -1.50;
	lboxx = 1.25;
	lboxy = 1.25;

%	Put a box around the flamelet
	plot([llcx, llcx+lboxx], [llcy, llcy], 'k-');
	plot([llcx, llcx+lboxx], [llcy+lboxy, llcy+lboxy], 'k-');
	plot([llcx+lboxx, llcx+lboxx], [llcy+lboxy, llcy], 'k-');
	plot([llcx, llcx], [llcy+lboxy, llcy], 'k-');

%	Plot the T-Z curve to represent the flamelet
	data = load('../data/lowStrain/lowStrain.CH4')';
	xnorm = data(:,1);
	ynorm = data(:,2);
	xnorm = xnorm - min(xnorm);
	ynorm = ynorm - min(ynorm);
	xnorm = xnorm./max(xnorm);
	ynorm = ynorm./max(ynorm);
	xnorm = xnorm*lboxx + llcx;
	ynorm  = ynorm*lboxy + llcy;
	plot(xnorm, ynorm);
	plot(xnorm(150), ynorm(150), 'o');

%	Formatting crap
	set(0, 'defaultTextInterpreter', 'LaTeX');

%	Label the flamelet
	text(llcx+1*lboxx/2, llcy+7*lboxy/8, '$\bar{p}_0$','FontSize', 14);
	text(llcx + lboxx/4, llcy+17*lboxy/16, '1D-Flame','FontSize', 14);
	text(llcx+3*lboxx/8, llcy-0.125, '$Z$', 'FontSize', 14);	
	text(llcx-lboxx/4, llcy+lboxy/2, '$T$', 'FontSize', 14);	
	text(llcx+7*lboxx/16, llcy+4*lboxy/8, '$\bar{Z}, \bar{T}_0$', 'FontSize', 14);

%	Add a second plot
	llcx2 = -3.;
	llcy2 =  0.25;
	lboxx2 = 1.25;
	lboxy2 = 1.25;
	
%	Put a box around the flamelet
	plot([llcx2, llcx2+lboxx2], [llcy2, llcy2], 'k-');
	plot([llcx2, llcx2+lboxx2], [llcy2+lboxy2, llcy2+lboxy2], 'k-');
	plot([llcx2+lboxx2, llcx2+lboxx2], [llcy2+lboxy2, llcy2], 'k-');
	plot([llcx2, llcx2], [llcy2+lboxy2, llcy2], 'k-');

%	Add a veritcal and horizontal line for overdoing it
	plot([llcx, llcx+lboxx], [ynorm(150), ynorm(150)], 'r:');
	plot([xnorm(150), xnorm(150)], [llcy, llcy+lboxy], 'r:');

%	Label the flamelet
	text(llcx2+3*lboxx2/8, llcy2+7*lboxy2/8, '$\bar{p}_0, \bar{T}_0$','FontSize', 14);
	text(llcx2+lboxx2/16, llcy2+17*lboxy2/16, 'Noble Gas-Air','FontSize', 14);
	text(llcx2+3*lboxx2/8, llcy2-0.125, '$Z$', 'FontSize', 14);	
	text(llcx2-lboxx2/4, llcy2+lboxy2/2, '$Y_{air}$', 'FontSize', 14);	
	text(llcx2+lboxx2, llcy2+lboxy2/2, '$Y_{noble}$', 'FontSize', 14);	
	text(llcx2+7*lboxx2/16, llcy2+2*lboxy2/8, '$\bar{Z}$', 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	PLOT THE INERT-AIR MIXTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	plot the mixing lines
	plot([llcx2, llcx2+lboxx2], [llcy2, llcy2+lboxy2], 'b-');
	plot([llcx2, llcx2+lboxx2], [llcy2+lboxy2, llcy2], 'm-');
	plot([llcx2+3*lboxx2/8], llcy2+3*lboxy2/8, 'ro');
	plot([llcx2+3*lboxx2/8], llcy2+5*lboxy2/8, 'ro');
	plot([llcx2+3*lboxx2/8, llcx2+3*lboxx2/8], [llcy2, llcy2 + lboxy2], 'r:');


%	Add a third plot to show where the mean state comes from
	llcxi = -1.25;
	llcyi =  -0.125;
	lboxxi = .25;
	lboxyi = .25;
	
%	Put a box around the flamelet
	plot([llcxi, llcxi+lboxxi], [llcyi, llcyi], 'k-');
	plot([llcxi, llcxi+lboxxi], [llcyi+lboxyi, llcyi+lboxyi], 'k-');
	plot([llcxi+lboxxi, llcxi+lboxxi], [llcyi+lboxyi, llcyi], 'k-');
	plot([llcxi, llcxi], [llcyi+lboxyi, llcyi], 'k-');

	plot([llcxi+lboxxi/2, llcx + lboxx], [llcyi+lboxyi/2, llcy+0*lboxy/2], 'k:');
	plot([llcxi+lboxxi/2, llcx + lboxx], [llcyi+lboxyi/2, llcy+lboxy], 'k:');

	plot([llcxi+lboxxi/2, llcx2 + lboxx2], [llcyi+lboxyi/2, llcy2+0*lboxy2/2], 'k:');
	plot([llcxi+lboxxi/2, llcx2 + lboxx2], [llcyi+lboxyi/2, llcy2+lboxy2], 'k:');

%	Formatting
	axis equal;
	axis([-3,3, -1.5, 1.5])
	axis off;

	print -depsc NozzleCartoon.eps
end
