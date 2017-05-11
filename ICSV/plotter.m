function[] = plotter()
	close all;
	load('processedFlucs.mat');

%	Flags
	raw_signals				= false;
	pressure_correction		= false;
	probe_cleaning			= false;
	final_spectra			= false;
	plot_characteristics	= false;
	plot_inlet_noise		= false;
	plot_outlet_chars		= false;
	plot_outlet_pi_spectra	= false;
	plot_outlet_pi_signal	= true;
	plot_outlet_p_spectra	= false;
	plot_inlet_signal		= false;
	plot_outlet_p_signal	= false;
	plot_SPL				= false;
	plot_good_SPL			= false;

%	Parameters
	L = 1;
	gamma = 1.4;
	gm1 = gamma - 1;
	gm1o2 = gm1/2;
	gogm1 = gamma/gm1;

%	Nozzle entrance condition
	Tbar_1 = mean(Tbar)
	pbar_1 = mean(pbar)
	rhobar_1 = mean(rhobar);
	Psibar_1 = -21;
	M_1 = 0.01;
	M_2 = 1.5;
	cbar_1 = sqrt(1.4*287*Tbar_1);
	c_p = 1004;

%	Stagnation condition
	Tbar_0 = Tbar_1/(1 + gm1o2*M_1*M_1)^(-1);
	pbar_0 = pbar_1/(1 + gm1o2*M_1*M_1)^(-gogm1);

%	Nozzle exit condition
	Tbar_2 = Tbar_0*(1 + gm1o2*M_2^2)^(-1);
	pbar_2 = pbar_0*(1 + gm1o2*M_2^2)^(-gogm1)

%	Plot the raw signals for evaluation
	Y = length(pprobe);
	t = linspace(0,N,N)'*dt;
	if (raw_signals)
		figure();
		subplot(2,2,1);
		plot(t, pprobe, 'LineWidth', 2);
		xlabel('time [ms]');
		ylabel('p_1 [Pa]');
		set(gca, 'FontSize', 16);
		subplot(2,2,2);
		plot(t, rhoprobe, 'LineWidth', 2);
		xlabel('time [ms]');
		ylabel('$\rho_1 [kg/m^3]$', 'Interpreter', 'LaTeX');
		set(gca, 'FontSize', 16);
		subplot(2,2,3);
		plot(t, vprobe, 'LineWidth', 2);
		xlabel('time [ms]');
		set(gca, 'FontSize', 16);
		ylabel('$w [m/s]$', 'Interpreter', 'LaTeX');
		subplot(2,2,4);
		plot(t, Zprobe, 'LIneWidth', 2);
		xlabel('time [ms]');
		ylabel('$Z [-]$', 'Interpreter', 'LaTeX');
		set(gca, 'FontSize', 16);
	end%(raw_signals)

%	Fit the pressure to a line, and then remove this linear trend from the signal
	pfit = polyfit(t, pprobe, 1);
	p_int = zeros(N,1);
	p_cor = zeros(N,1);
	for i = 1:N
		p_int(i) = polyval(pfit, t(i));
		p_cor(i) = pprobe(i) - p_int(i);
	end
	p_cor = p_cor - mean(p_cor);

%	Make two plots to demonstrate the pressure correction procedure
	if (pressure_correction)
		figure();
		plot(t, pprobe, 'b-', 'LineWidth', 2);
		hold on;
		plot(t, p_int, 'r--', 'LineWidth', 6);
		xlabel('time [ms]');
		ylabel('p [Pa]');
		aa = legend('Raw signal', 'Linear trend');
		set(aa, 'FontSize', 14);
		set(gca, 'FontSize', 14);

		figure();
		plot(t, p_cor, 'b-', 'LineWidth', 2);
		hold on;
		xlabel('time [ms]');
		ylabel('p [Pa]');
		set(gca, 'FontSize', 14);
	end
	
	[pprobe] = cleanProbes(t, p_cor, probe_cleaning);
%	[sigma_probe] = cleanProbes(t, sigma_probe);
	[wprobe] = cleanProbes(t, wprobe, probe_cleaning);
	[rhoprobe] = cleanProbes(t, rhoprobe, probe_cleaning);
	[Zprobe, t_smpl] = cleanProbes(t, Zprobe, probe_cleaning);
	[N, Nprobes] = size(pprobe);
	dt = 250*5E-8;
	omega = 1/(N*dt)*[-N/2:N/2-1];

%	Build raw probe inputs
%	Construct the entropy fluctuations
	sigma_probe = pprobe/(gamma*pbar_1) - rhoprobe/rhobar_1 - Psibar_1*Zprobe;
%	Window the probes so that they look nice

	p_DM = zeros(N,Nprobes);
	p_MC = zeros(N,Nprobes);
	s_DM = zeros(N,Nprobes);
	s_MC = zeros(N,Nprobes);
	z_DM = zeros(N,Nprobes);
	z_MC = zeros(N,Nprobes);

	if (plot_characteristics)
		pstar_plot 			= figure();
		set(pstar_plot, 'Position', [0 0 1500, 250]);
		wstar_plot 		= figure();
		set(wstar_plot, 'Position', [0 0 1500, 250]);
		rhostar_plot	= figure();
		set(rhostar_plot, 'Position', [0 0 1500, 250]);
		Pseye_plot		= figure();
		set(Pseye_plot, 'Position', [0 0 1500, 250]);
		pi_plot 		= figure();
		set(pi_plot, 'Position', [0 0 1500, 250]);
		sigma_plot		= figure();
		set(sigma_plot, 'Position', [0 0 1500, 250]);
		xi_plot			= figure();
		set(xi_plot, 'Position', [0 0 1500, 250]);
	end%plot_characteristics

	for probe = 1:Nprobes
		p1_star(:,probe) 	= pprobe(:,probe)/(gamma*pbar_1);
		w1_star(:,probe) 	= wprobe(:,probe)/cbar_1;
		rho1_star(:,probe)	= rhoprobe(:,probe)/rhobar_1;
		Pseye1(:,probe)		= Zprobe(:,probe).*Psibar_1;

		pi_plus1(:,probe) 	= 0.5*(p1_star(:,probe) + w1_star(:,probe));
		sigma1(:,probe)		= p1_star(:,probe) - rho1_star(:,probe) - Pseye1(:,probe);
		xi1(:,probe)		= Zprobe(:,probe);

		if (plot_characteristics)
			figure(pstar_plot);
			subplot(1,7,probe);
			plot(t_smpl(:,probe),p1_star(:,probe), 'LineWidth', 2);
			xlabel('Time [ms]');
			if (probe == 1)
				ylabel('$p^* \slash \gamma \bar{p}_a$', 'Interpreter', 'LaTeX');
			end
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 14);

			figure(wstar_plot);
			subplot(1,7,probe);
			plot(t_smpl(:,probe),w1_star(:,probe), 'LineWidth', 2);
			xlabel('Time [ms]');
			if (probe == 1)
				ylabel('$w^* \slash \bar{c}_a$', 'Interpreter', 'LaTeX');
			end
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 14);

			figure(pi_plot);
			subplot(1,7,probe);
			plot(t_smpl(:,probe),pi_plus1(:,probe), 'LineWidth', 2);
			xlabel('Time [ms]');
			if (probe == 1)
				ylabel('$\pi^+_a$', 'Interpreter', 'LaTeX');
			end
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 14);

			figure(rhostar_plot);
			subplot(1,7,probe);
			plot(t_smpl(:,probe),rho1_star(:,probe), 'LineWidth', 2);
			xlabel('Time [ms]');
			if (probe == 1)
				ylabel('$\rho^* \slash \bar{\rho}_1$', 'Interpreter', 'LaTeX');
			end
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 14);

			figure(Pseye_plot);
			subplot(1,7,probe);
			plot(t_smpl(:,probe),Pseye1(:,probe), 'LineWidth', 2);
			xlabel('Time [ms]');
			if (probe == 1)
				ylabel('$\bar{\Psi}_a \xi^*$', 'Interpreter', 'LaTeX');
			end
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 14);

			figure(sigma_plot);
			subplot(1,7,probe);
			plot(t_smpl(:,probe),sigma1(:,probe), 'LineWidth', 2);
			xlabel('Time [ms]');
			if (probe == 1)
				ylabel('$\sigma_a$', 'Interpreter', 'LaTeX');
			end
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 14);

			figure(xi_plot);
			subplot(1,7,probe);
			plot(t_smpl(:,probe),xi1(:,probe), 'LineWidth', 2);
			xlabel('Time [ms]');
			if (probe == 1)
				ylabel('$\xi_a$', 'Interpreter', 'LaTeX');
			end
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 14);
		end%plot_characteristics
	end%probe = 1:Nprobes
	if (plot_characteristics)
		tightfig(pstar_plot);
		tightfig(wstar_plot);
		tightfig(rhostar_plot);
		tightfig(Pseye_plot);
		tightfig(xi_plot);
		tightfig(sigma_plot);
		tightfig(pi_plot);
	end%(plot_characteristics)

	if (plot_inlet_noise)
		a1 = figure();
		set(a1, 'Position', [ 0 0 1500 250]);	
		a2 = figure();
		set(a2, 'Position', [ 0 0 1500 250]);	
		a3 = figure();
		set(a3, 'Position', [ 0 0 1500 250]);	
	end

	N = length(pprobe(:,1));
	omega_smpl = 1/(N*dt)*[-N/2:N/2-1];
	He = omega_smpl*L/cbar_1;;
	for probe = 1:Nprobes
		pprobe_hat(:,probe) = fftshift(fft(pprobe(:,probe)));
		pprobe_hat(N/2+1,probe) = NaN;
		SPL_1_hat(:,probe) = 10*log10((pprobe_hat(:,probe)/20E-6).*conj(pprobe_hat(:,probe)/20E-6))/sqrt(2);
		SPL_1_hat(N/2+1,probe) = NaN;
		if (plot_inlet_noise)
			figure(a1);
			subplot(1,8,probe);
			plot(t_smpl(:,probe), pprobe(:,probe),'LineWidth', 2);
			title(strcat('Sample #', num2str(probe)));
			xlabel('t [ms]');
			set(gca, 'FontSize', 16);
			if (probe == 1)
				ylabel('$p$ [Pa]', 'Interpreter', 'LaTeX');
			end
			xlim([min(t_smpl(:,probe)) max(t_smpl(:,probe))]);

			figure(a2);
			subplot(1,8,probe);
			semilogy(omega_smpl, abs(pprobe_hat(:,probe)), 'LineWidth', 2);
			title(strcat('Sample #', num2str(probe)));
			xlabel('f [Hz]');
			xlim([-2000,2000]);
			ylim([1E3, 1E6]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				ylabel('$\hat{p}$ [Pa-Hz]', 'Interpreter', 'LaTeX');
			end

			figure(a3);
			subplot(1,8,probe);
			plot(omega_smpl, SPL_1_hat(:,probe), 'LineWidth', 2);
			title(strcat('Sample #', num2str(probe)));
			set(gca, 'FontSize', 16);
			xlabel('f [Hz]');
			ylim([110,150]);
			xlim([-2000,2000]);
			if (probe == 1)
				ylabel('SPL [dB]', 'Interpreter', 'LaTeX');
			end
		end%(plot_inlet_noise)
	end

	pprobe_hat_bar = mean(abs(pprobe_hat)');
	SPL_1_hat_bar = mean(SPL_1_hat')';

	if (plot_inlet_noise)
		figure(a1);
		set(gca ,'FontSize', 16);
		tightfig();

		figure(a2);
		subplot(1,8,8);
		semilogy(omega_smpl, pprobe_hat_bar, 'LineWidth', 2);
		xlabel('f [Hz]');
		xlim([-2000,2000]);
		ylim([1E3, 1E6]);
		set(gca ,'FontSize', 16);
		title('Average Respone');
		tightfig();

		figure(a3);
		subplot(1,8,8);
		plot(omega_smpl, SPL_1_hat_bar, 'LineWidth', 2);
		xlabel('f [Hz]');
		xlim([-2000,2000]);
		set(gca, 'FontSize', 16);
		title('Average Response');
		tightfig();
	end%plot_inlet_noise

%	Loop over all the probes and apply the model
	for probe = 1:Nprobes
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%				DIRECT NOISE				  %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		pi_plus_hat_1(:,probe) = fftshift(fft(pi_plus1(:,probe)));

%		Load the data from the Helmholtz dependent solver
		load('Acoustic.mat')									% Load data
		OMEGA = OMEGA';											% Force column vectors
		modulusTransferFn = abs(TRANS(:,1));					% Peel off the downstream acoustic transfer function from the TRANS array
		minusmodulusTransferFn = abs(TRANS(:,2));
		modulus = spline(OMEGA, modulusTransferFn);				% Spline the transfer function modulus to allow for continuous queries
		minusmodulus = spline(OMEGA, minusmodulusTransferFn);
		phaseTransferFn = PHASE(:,1);							% Peel off the downstream acoustic transfer function from the PHASE array
		minusphaseTransferFn = PHASE(:,2);
		phase 	= spline(OMEGA, phaseTransferFn);
		minusphase	= spline(OMEGA, minusphaseTransferFn);

%		Allocate some arrays
		pi_plus_hat_2_DM	= zeros(N,Nprobes);						% Store the dimensionless downstream pressure wave computed with Duran-Moreau style model
		pi_plus_hat_2_MC	= zeros(N,Nprobes);
		p_2_hat_amp_DM 			= zeros(N,Nprobes);						% Store the dimensional downstream pressure wave computed with Duran-Moreau style model
		p_2_hat_amp_MC 			= zeros(N,Nprobes);						% Store the dimensional downstream pressure wave computed with Marble-Candel style model

%		Apply the transfer functions to the incoming data			
		for i = 1:N
			if (abs(He(i)) < 2)													% Only act on waves w He < 2
				p_phase_1(i,probe)			= atan2(imag(pi_plus_hat_1(i,probe)), real(pi_plus_hat_1(i,probe)));
				pi_plus_hat_2_DM(i,probe)	= abs(pi_plus_hat_1(i,probe))*ppval(modulus, abs(He(i)));
				pi_plus_hat_2_MC(i,probe)	= abs(pi_plus_hat_1(i,probe))*ppval(modulus, 0);
				pi_minus_hat_2_DM(i,probe)	= abs(pi_plus_hat_1(i,probe))*ppval(minusmodulus, abs(He(i)));
				pi_minus_hat_2_MC(i,probe)	= abs(pi_plus_hat_1(i,probe))*ppval(minusmodulus, 0);
				pi_total_hat_DM(i,probe)	= pi_plus_hat_2_DM(i,probe)*(cos(pi_phase_2_DM(i,probe)) + sqrt(-1)*pi_phase_2_DM(i,probe));
				pi_total_hat_MC(i,probe)	= pi_plus_hat_2_MC(i,probe)*(cos(pi_phase_2_MC(i,probe)) + sqrt(-1)*pi_phase_2_MC(i,probe));
				p_2_hat_amp_DM(i,probe) 	= abs(pi_plus_hat_2_DM(i,probe))*gamma*pbar_2;			% Scale the dimensionless output to make it dimensional MC		
				p_2_hat_amp_MC(i,probe) 	= abs(pi_plus_hat_2_MC(i,probe))*gamma*pbar_2;			% Scale the dimensionless output to make it dimensional DM	
				p_phase_2_DM(i,probe)		= p_phase_1(i,probe) + ppval(phase, abs(He(i)));
				p_phase_2_MC(i,probe)		= p_phase_1(i,probe) + ppval(phase, 0);
%				p_phase_2_MC(i,probe) = 0;
			else
				p_phase_1(i,probe) = 0;
				p_phase_2_DM(i,probe) = 0;
				p_phase_2_MC(i,probe) = 0;
			end
		end

%		Compute SPL from pressure fluctuation levels	
%		SPL_2_hat_DM(:,probe) 		= 10*log10((p_2_hat_amp_DM(:,probe)/20E-6).*conj(p_2_hat_amp_DM(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to DM pressure
%		SPL_2_hat_MC(:,probe)		= 10*log10((p_2_hat_amp_MC(:,probe)/20E-6).*conj(p_2_hat_amp_MC(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to MC pressure
%		SPL_2_hat_DM(N/2+1,probe)	= NaN;																					% Null out the mean for aesthetics DM
%		SPL_2_hat_MC(N/2+1,probe) 	= NaN;																					% Null out the mean for aesthetics MC

%		Build some arrays for output & plotting
%		p_pi_hat_DM(:,probe) 		= pi_plus_hat_2_DM(:,probe);
%		p_pi_hat_MC(:,probe)		= pi_plus_hat_2_MC(:,probe);
		p_pi_hat_DM(:,probe)		= abs(pi_plus_hat_2_DM(:,probe)).*(cos(p_phase_2_DM(:,probe)) + sqrt(-1).*sin(p_phase_2_DM(:,probe)));
		p_pi_hat_MC(:,probe)		= abs(pi_plus_hat_2_MC(:,probe)).*(cos(p_phase_2_MC(:,probe)) + sqrt(-1)*sin(p_phase_2_MC(:,probe)));
		p_pi_DM(:,probe)			= ifft(fftshift(p_pi_hat_DM(:,probe)));
		p_pi_MC(:,probe)			= ifft(fftshift(p_pi_hat_MC(:,probe)));
%		p_p_hat_DM(:,probe)			= p_2_hat_DM(:,probe);
%		p_p_hat_MC(:,probe)			= p_2_hat_MC(:,probe);
		p_p_hat_DM(:,probe)			= p_2_hat_amp_DM(:,probe).*(cos(p_phase_2_DM(:,probe)) + sqrt(-1)*sin(p_phase_2_DM(:,probe)));
		p_p_hat_MC(:,probe)			= p_2_hat_amp_MC(:,probe).*(cos(p_phase_2_MC(:,probe)) + sqrt(-1)*sin(p_phase_2_MC(:,probe)));
%%%%		
		SPL_2_hat_DM(:,probe) 	= 10*log10((p_p_hat_DM(:,probe)/20E-6).*conj(p_p_hat_DM(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to DM pressure
		SPL_2_hat_MC(:,probe)	= 10*log10((p_p_hat_MC(:,probe)/20E-6).*conj(p_p_hat_MC(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to MC pressure
		SPL_2_hat_DM(N/2+1,probe)	= NaN;														% Null out the mean for aesthetics DM
		SPL_2_hat_MC(N/2+1,probe) = NaN;														% Null out the mean for aesthetics MC
%%%%
		p_p_DM(:,probe)			= ifft(fftshift(p_p_hat_DM(:,probe)));
		p_p_MC(:,probe)			= ifft(fftshift(p_p_hat_MC(:,probe)));
		p_SPL_DM(:,probe) 	= SPL_2_hat_DM(:,probe);
		p_SPL_MC(:,probe) 	= SPL_2_hat_MC(:,probe);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%				ENTROPY NOISE				  %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		sigma_hat_1(:,probe) = fftshift(fft(sigma1(:,probe)));

%		Load the data from the Helmholtz dependent solver
		load('Entropy.mat')										% Load data
		OMEGA = OMEGA';											% Force column vectors
		modulusTransferFn = abs(TRANS(:,1));					% Peel off the downstream acoustic transfer function from the TRANS array
		phaseTransferFn = PHASE(:,1);							% Peel off the downstream acoustic transfer function from the PHASE array
		modulus = spline(OMEGA, modulusTransferFn);				% Spline the transfer function modulus to allow for continuous queries
		phase 	= spline(OMEGA, phaseTransferFn);

%		Allocate some arrays
		pi_plus_hat_2_DM	= zeros(N,Nprobes);						% Store the dimensionless downstream pressure wave computed with Duran-Moreau style model
		pi_plus_hat_2_MC	= zeros(N,Nprobes);						% Store the dimensionless downstream pressure wave computed with Marble-Candel style model
		p_2_hat_amp_DM 		= zeros(N,Nprobes);						% Store the dimensional downstream pressure wave computed with Duran-Moreau style model
		p_2_hat_amp_MC 		= zeros(N,Nprobes);						% Store the dimensional downstream pressure wave computed with Marble-Candel style model

%		Apply the transfer functions to the incoming data			
		for i = 1:N
			if (abs(He(i)) < 2)																	% Only act on waves w He < 2
				s_phase_1(i,probe)			= atan2(imag(sigma_hat_1(i,probe)), real(sigma_hat_1(i,probe)));
				pi_plus_hat_2_DM(i,probe) 	= abs(sigma_hat_1(i,probe))*ppval(modulus, abs(He(i)));	% Stretch the modulus by the appropriate factor for current He
				pi_plus_hat_2_MC(i,probe)	= abs(sigma_hat_1(i,probe))*ppval(modulus, 0);			% Stretch the modulus by the factor for He = 0
				p_2_hat_amp_DM(i,probe) 	= abs(pi_plus_hat_2_DM(i,probe))*gamma*pbar_2;			% Scale the dimensionless output to make it dimensional MC		
				p_2_hat_amp_MC(i,probe) 	= abs(pi_plus_hat_2_MC(i,probe))*gamma*pbar_2;			% Scale the dimensionless output to make it dimensional DM	
				s_phase_2_DM(i,probe)		= s_phase_1(i,probe) + ppval(phase, abs(He(i)));
				s_phase_2_MC(i,probe)		= s_phase_1(i,probe) + ppval(phase, 0);
%				s_phase_2_MC(i,probe) = 0;
			else
				s_phase_1(i,probe) = 0;
				s_phase_2_DM(i,probe) = 0;
				s_phase_2_MC(i,probe) = 0;
			end
		end

%		Compute SPL from pressure fluctuation levels	
%		SPL_2_hat_DM(:,probe)	= 10*log10(p_2_hat_amp_DM(:,probe)/20E-6.*conj(p_2_hat_amp_DM(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to DM pressure
%		SPL_2_hat_MC(:,probe)	= 10*log10(p_2_hat_amp_MC(:,probe)/20E-6.*conj(p_2_hat_amp_MC(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to MC pressure
%		SPL_2_hat_DM(N/2+1,probe)	= NaN;																					% Null out the mean for aesthetics DM
%		SPL_2_hat_MC(N/2+1,probe) 	= NaN;																					% Null out the mean for aesthetics MC

%		Build some arrays for output & plotting
		s_pi_hat_DM(:,probe)		= abs(pi_plus_hat_2_DM(:,probe)).*(cos(s_phase_2_DM(:,probe)) + sqrt(-1).*sin(s_phase_2_DM(:,probe)));
		s_pi_hat_MC(:,probe)		= abs(pi_plus_hat_2_MC(:,probe)).*(cos(s_phase_2_MC(:,probe)) + sqrt(-1)*sin(s_phase_2_MC(:,probe)));
%		s_pi_hat_DM(:,probe) 		= pi_plus_hat_2_DM(:,probe);
%		s_pi_hat_MC(:,probe)		= pi_plus_hat_2_MC(:,probe);
		s_pi_DM(:,probe)			= ifft(fftshift(s_pi_hat_DM(:,probe)));
		s_pi_MC(:,probe)			= ifft(fftshift(s_pi_hat_MC(:,probe)));
		s_p_hat_DM(:,probe)			= p_2_hat_amp_DM(:,probe).*(cos(s_phase_2_DM(:,probe)) + sqrt(-1)*sin(s_phase_2_DM(:,probe)));
		s_p_hat_MC(:,probe)			= p_2_hat_amp_MC(:,probe).*(cos(s_phase_2_MC(:,probe)) + sqrt(-1)*sin(s_phase_2_MC(:,probe)));
%		s_p_hat_DM(:,probe)			= p_2_hat_amp_DM(:,probe);
%		s_p_hat_MC(:,probe)			= p_2_hat_amp_MC(:,probe);
%%%%		
		SPL_2_hat_DM(:,probe) 	= 10*log10((s_p_hat_DM(:,probe)/20E-6).*conj(s_p_hat_DM(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to DM pressure
		SPL_2_hat_MC(:,probe)	= 10*log10((s_p_hat_MC(:,probe)/20E-6).*conj(s_p_hat_MC(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to MC pressure
		SPL_2_hat_DM(N/2+1,probe)	= NaN;														% Null out the mean for aesthetics DM
		SPL_2_hat_MC(N/2+1,probe) = NaN;														% Null out the mean for aesthetics MC
%%%%
		s_p_DM(:,probe)				= ifft(fftshift(s_p_hat_DM(:,probe)));
		s_p_MC(:,probe)				= ifft(fftshift(s_p_hat_MC(:,probe)));
		s_SPL_DM(:,probe) 			= SPL_2_hat_DM(:,probe);
		s_SPL_MC(:,probe) 			= SPL_2_hat_MC(:,probe);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%			COMPOSITION NOISE				  %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		xi_hat_1(:,probe) = fftshift(fft(xi1(:,probe)));
%		Load the data from the Helmholtz dependent solver
		load('Composition.mat')									% Load data
		OMEGA = OMEGA';											% Force column vectors
		modulusTransferFn = abs(TRANS(:,1));					% Peel off the downstream acoustic transfer function from the TRANS array
		phaseTransferFn = PHASE(:,1);							% Peel off the downstream acoustic transfer function from the PHASE array
		modulus = spline(OMEGA, modulusTransferFn);				% Spline the transfer function modulus to allow for continuous queries
		phase 	= spline(OMEGA, phaseTransferFn);
	
		if (1 == 2)
		if (probe == 1)
			figure();
			plot(OMEGA, phaseTransferFn*180/pi);
			end
		end

%		Allocate some arrays
		p_phase_1			= zeros(N,Nprobes);
		z_phase_2_DM		= zeros(N,Nprobes);
		z_phase_2_MC		= zeros(N,Nprobes);
		pi_plus_hat_2_DM	= zeros(N,Nprobes);						% Store the dimensionless downstream pressure wave computed with Duran-Moreau style model
		pi_plus_hat_2_MC	= zeros(N,Nprobes);						% Store the dimensionless downstream pressure wave computed with Marble-Candel style model
		p_2_hat_amp_DM 		= zeros(N,Nprobes);						% Store the dimensional downstream pressure wave computed with Duran-Moreau style model
		p_2_hat_amp_MC 		= zeros(N,Nprobes);						% Store the dimensional downstream pressure wave computed with Marble-Candel style model

%		Apply the transfer functions to the incoming data			
		for i = 1:N
			if (abs(He(i)) < 2)																	% Only act on waves w He < 2
				z_phase_1(i,probe)			= atan2(imag(xi_hat_1(i,probe)), real(xi_hat_1(i,probe)));
				pi_plus_hat_2_DM(i,probe) 	= abs(xi_hat_1(i,probe))*ppval(modulus, abs(He(i)));		% Stretch the modulus by the appropriate factor for current He
				pi_plus_hat_2_MC(i,probe)	= abs(xi_hat_1(i,probe))*ppval(modulus, 0);				% Stretch the modulus by the factor for He = 0
				p_2_hat_amp_DM(i,probe) 	= abs(pi_plus_hat_2_DM(i,probe))*gamma*pbar_2;			% Scale the dimensionless output to make it dimensional MC		
				p_2_hat_amp_MC(i,probe) 	= abs(pi_plus_hat_2_MC(i,probe))*gamma*pbar_2;			% Scale the dimensionless output to make it dimensional DM	
				z_phase_2_DM(i,probe)		= z_phase_1(i,probe) + ppval(phase, abs(He(i)));
				z_phase_2_MC(i,probe)		= z_phase_1(i,probe) + ppval(phase, 0);
%				z_phase_2_MC(i,probe) = 0;
			else
				z_phase_1(i,probe) = 0;
				z_phase_2_DM(i,probe) = 0;
				z_phase_2_MC(i,probe) = 0;
			end
		end
		if (1 == 2)
			if (probe == 1)
			hold on;
			plot(He, 180/pi*(z_phase_2_DM(:,1) - z_phase_1(:,1)), '-o', 'LineWidth', 2);
			xlabel('He');
			ylabel('phase lag');
			xlim([-2,2]);
			end
		end

%		Compute SPL from pressure fluctuation levels	
%		SPL_2_hat_DM(:,probe) 	= 10*log10((p_2_hat_amp_DM(:,probe)/20E-6).*conj(p_2_hat_amp_DM(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to DM pressure
%		SPL_2_hat_MC(:,probe)	= 10*log10((p_2_hat_amp_MC(:,probe)/20E-6).*conj(p_2_hat_amp_MC(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to MC pressure
%		SPL_2_hat_DM(N/2+1,probe)	= NaN;														% Null out the mean for aesthetics DM
%		SPL_2_hat_MC(N/2+1,probe) = NaN;														% Null out the mean for aesthetics MC

%		Build some arrays for output & plotting
%		z_pi_hat_DM(:,probe) 		= pi_plus_hat_2_DM(:,probe);
%		z_pi_hat_MC(:,probe)		= pi_plus_hat_2_MC(:,probe);
		z_pi_hat_DM(:,probe)		= abs(pi_plus_hat_2_DM(:,probe)).*(cos(z_phase_2_DM(:,probe)) + sqrt(-1).*sin(z_phase_2_DM(:,probe)));
		z_pi_hat_MC(:,probe)		= abs(pi_plus_hat_2_MC(:,probe)).*(cos(z_phase_2_MC(:,probe)) + sqrt(-1)*sin(z_phase_2_MC(:,probe)));
		z_pi_DM(:,probe)			= ifft(fftshift(z_pi_hat_DM(:,probe)));
		z_pi_MC(:,probe)			= ifft(fftshift(z_pi_hat_MC(:,probe)));
%		z_p_hat_DM(:,probe)			= p_2_hat_amp_DM(:,probe);
%		z_p_hat_MC(:,probe)			= p_2_hat_amp_MC(:,probe);
		z_p_hat_DM(:,probe)			= p_2_hat_amp_DM(:,probe).*(cos(z_phase_2_DM(:,probe)) + sqrt(-1)*sin(z_phase_2_DM(:,probe)));
		z_p_hat_MC(:,probe)			= p_2_hat_amp_MC(:,probe).*(cos(z_phase_2_MC(:,probe)) + sqrt(-1)*sin(z_phase_2_MC(:,probe)));

%%%%		
		SPL_2_hat_DM(:,probe) 	= 10*log10((z_p_hat_DM(:,probe)/20E-6).*conj(z_p_hat_DM(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to DM pressure
		SPL_2_hat_MC(:,probe)	= 10*log10((z_p_hat_MC(:,probe)/20E-6).*conj(z_p_hat_MC(:,probe)/20E-6))/sqrt(2);		% Apply the SPL definition to MC pressure
		SPL_2_hat_DM(N/2+1,probe)	= NaN;														% Null out the mean for aesthetics DM
		SPL_2_hat_MC(N/2+1,probe) = NaN;														% Null out the mean for aesthetics MC
%%%%
		
		z_p_DM(:,probe)				= ifft(fftshift(z_p_hat_DM(:,probe)));
		z_p_MC(:,probe)				= ifft(fftshift(z_p_hat_MC(:,probe)));
		z_SPL_DM(:,probe) 			= SPL_2_hat_DM(:,probe);
		z_SPL_MC(:,probe) 			= SPL_2_hat_MC(:,probe);
	end%for probe

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%				TOTAL NOISE					  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	t_p_hat_DM = p_p_hat_DM + s_p_hat_DM + z_p_hat_DM;
	t_p_hat_MC = p_p_hat_MC + s_p_hat_MC + z_p_hat_MC;
%	t_p_hat_DM = s_p_hat_DM + z_p_hat_DM;
%	t_p_hat_MC = s_p_hat_MC + z_p_hat_MC;

	for probe = 1:Nprobes
		SPL_2_hat_DM(:,probe)	= 10*log10((t_p_hat_DM(:,probe)/20E-6).*conj(t_p_hat_DM(:,probe)/20E-6))/sqrt(2);
		SPL_2_hat_MC(:,probe)	= 10*log10((t_p_hat_MC(:,probe)/20E-6).*conj(t_p_hat_MC(:,probe)/20E-6))/sqrt(2);
		SPL_2_hat_DM(N/2+1,probe)	= NaN;
		SPL_2_hat_MC(N/2+1,probe)	= NaN;

		t_SPL_DM(:,probe)			= SPL_2_hat_DM(:,probe);
		t_SPL_MC(:,probe)			= SPL_2_hat_MC(:,probe);
	end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%				POST-PROCESSING				  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	Do averaging for post-processing	
	p_SPL_bar_DM 		= mean(p_SPL_DM')';
	p_SPL_bar_MC 		= mean(p_SPL_MC')';
	s_SPL_bar_DM 		= mean(s_SPL_DM')';
	s_SPL_bar_MC 		= mean(s_SPL_MC')';
	z_SPL_bar_DM 		= mean(z_SPL_DM')';
	z_SPL_bar_MC 		= mean(z_SPL_MC')';
	t_SPL_bar_DM 		= mean(t_SPL_DM')';
	t_SPL_bar_MC 		= mean(t_SPL_MC')';

	if (plot_outlet_pi_spectra)
		c1 = figure();
		set(c1, 'Position', [0 0 1500 250]);
		c2 = figure();
		set(c2, 'Position', [0 0 1500 250]);
		c3 = figure();
		set(c3, 'Position', [0 0 1500 250]);
		for probe = 1:Nprobes
			figure(c1);
			subplot(1,8,probe);
			semilogy(omega_smpl, abs(p_pi_hat_MC(:,probe)), 'k-', 'LineWidth', 2);
			hold on;
			semilogy(omega_smpl, abs(p_pi_hat_DM(:,probe)), 'b-', 'LineWidth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('$\hat{\pi}^+_b$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([1E1, 1E7]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end

			figure(c2);
			subplot(1,8,probe);
			semilogy(omega_smpl, abs(s_pi_hat_MC(:,probe)), 'k-', 'LineWidth', 2);
			hold on;
			semilogy(omega_smpl, abs(s_pi_hat_DM(:,probe)), 'b-', 'LineWidth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('$\hat{\pi}^+_b$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([1E1, 1E7]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end

			figure(c3);
			subplot(1,8,probe);
			semilogy(omega_smpl, abs(z_pi_hat_MC(:,probe)), 'k-', 'lineWidth', 2);
			hold on;
			semilogy(omega_smpl, abs(z_pi_hat_DM(:,probe)), 'b-', 'LineWIdth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('$\hat{\pi}^+_b$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([1E1, 1E7]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end
		end%for probe

		figure(c1);
		subplot(1,8,8);
		semilogy(omega_smpl, mean(abs(p_pi_hat_MC)')', 'k-', 'LineWidth', 2);
		hold on;
		semilogy(omega_smpl, mean(abs(p_pi_hat_DM)')', 'b-', 'LineWidth', 2);
		xlabel('f [Hz]');
		title('Average');
		xlim([-2000, 2000]);
		ylim([1E1, 1E7]);
		set(gca, 'FontSize', 16);
		tightfig();

		figure(c2);
		subplot(1,8,8);
		semilogy(omega_smpl, mean(abs(s_pi_hat_MC)')', 'k-', 'LineWidth', 2);
		hold on;
		semilogy(omega_smpl, mean(abs(s_pi_hat_DM)')', 'b-', 'LineWidth', 2);
		title('Average');
		xlabel('f [Hz]');
		xlim([-2000, 2000]);
		ylim([1E1, 1E7]);
		set(gca, 'FontSize', 16);
		tightfig();

		figure(c3);
		subplot(1,8,8);
		semilogy(omega_smpl, mean(abs(z_pi_hat_MC)')', 'k-', 'LineWidth', 2);
		hold on;
		semilogy(omega_smpl, mean(abs(z_pi_hat_DM)')', 'b-', 'LineWidth', 2);
		title('Average');
		xlabel('f [Hz]');
		xlim([-2000, 2000]);
		ylim([1E1, 1E7]);
		set(gca, 'FontSize', 16);
		tightfig();


	end%(plot_outlet_spectra)

	if (plot_outlet_pi_signal)
		d1 = figure();
		set(d1, 'Position', [0 0 1500 250]);
		d2 = figure();
		set(d2, 'Position', [0 0 1500 250]);
		d3 = figure();
		set(d3, 'Position', [0 0 1500 250]);
		for probe = 1:Nprobes
			figure(d1);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), p_pi_MC(:,probe), 'k-', 'LineWidth', 2);
			hold on;
			plot(t_smpl(:,probe), p_pi_DM(:,probe), 'b-', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$\pi^+_b$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end

			figure(d2);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), s_pi_MC(:,probe), 'k-', 'LineWidth', 2);
			hold on;
			plot(t_smpl(:,probe), s_pi_DM(:,probe), 'b-', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$\pi^+_b$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end
			
			figure(d3);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), z_pi_MC(:,probe), 'k-', 'LineWidth', 2);
			hold on;
			plot(t_smpl(:,probe), z_pi_DM(:,probe), 'b-', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$\pi^+_b$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end
		end%for probe
	end%(plot_outlet_spectra)



	if (plot_outlet_p_spectra)
		e1 = figure();
		set(e1, 'Position', [0 0 1500 250]);
		e2 = figure();
		set(e2, 'Position', [0 0 1500 250]);
		e3 = figure();
		set(e3, 'Position', [0 0 1500 250]);
		for probe = 1:Nprobes
			figure(e1);
			subplot(1,8,probe);
			semilogy(omega_smpl, abs(p_p_hat_MC(:,probe)), 'k-', 'LineWidth', 2);
			hold on;
			semilogy(omega_smpl, abs(p_p_hat_DM(:,probe)), 'b-', 'LineWidth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('$\hat{p}$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([1E2, 1E7]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end

			figure(e2);
			subplot(1,8,probe);
			semilogy(omega_smpl, abs(s_p_hat_MC(:,probe)), 'k-', 'LineWidth', 2);
			hold on;
			semilogy(omega_smpl, abs(s_p_hat_DM(:,probe)), 'b-', 'LineWidth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('$\hat{p}$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([1E2, 1E7]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end

			figure(e3);
			subplot(1,8,probe);
			semilogy(omega_smpl, abs(z_p_hat_MC(:,probe)), 'k-', 'lineWidth', 2);
			hold on;
			semilogy(omega_smpl, abs(z_p_hat_DM(:,probe)), 'b-', 'LineWIdth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('$\hat{p}$', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([1E2, 1E7]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end
		end%for probe

		figure(e1);
		subplot(1,8,8);
		semilogy(omega_smpl, mean(abs(p_p_hat_MC)')', 'k-', 'LineWidth', 2);
		hold on;
		semilogy(omega_smpl, mean(abs(p_p_hat_DM)')', 'b-', 'LineWidth', 2);
		xlabel('f [Hz]');
		title('Average');
		xlim([-2000, 2000]);
		ylim([1E2, 1E7]);
		set(gca, 'FontSize', 16);
		tightfig();

		figure(e2);
		subplot(1,8,8);
		semilogy(omega_smpl, mean(abs(s_p_hat_MC)')', 'k-', 'LineWidth', 2);
		hold on;
		semilogy(omega_smpl, mean(abs(s_p_hat_DM)')', 'b-', 'LineWidth', 2);
		title('Average');
		xlabel('f [Hz]');
		xlim([-2000, 2000]);
		ylim([1E2, 1E7]);
		set(gca, 'FontSize', 16);
		tightfig();

		figure(e3);
		subplot(1,8,8);
		semilogy(omega_smpl, mean(abs(z_p_hat_MC)')', 'k-', 'LineWidth', 2);
		hold on;
		semilogy(omega_smpl, mean(abs(z_p_hat_DM)')', 'b-', 'LineWidth', 2);
		title('Average');
		xlabel('f [Hz]');
		xlim([-2000, 2000]);
		ylim([1E2, 1E7]);
		set(gca, 'FontSize', 16);
		tightfig();


	end%(plot_outlet_spectra)

	if (plot_outlet_p_signal)
		f1 = figure();
		set(f1, 'Position', [0 0 1500 250]);
		f2 = figure();
		set(f2, 'Position', [0 0 1500 250]);
		f3 = figure();
		set(f3, 'Position', [0 0 1500 250]);
		for probe = 1:Nprobes
			figure(f1);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), real(p_p_MC(:,probe)), 'k-', 'LineWidth', 2);
			hold on;
			plot(t_smpl(:,probe), real(p_p_DM(:,probe)), 'b:', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$p$ [Pa]', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
%			ylim([-10000,10000]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'NorthEast');
			end

			figure(f2);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), real(s_p_MC(:,probe)), 'k-', 'LineWidth', 2);
			hold on;
			plot(t_smpl(:,probe), real(s_p_DM(:,probe)), 'b:', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$p$ [Pa]', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			ylim([-10000,10000]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'NorthEast');
			end
			
			figure(f3);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), real(z_p_MC(:,probe)), 'k-', 'LineWidth', 2);
			hold on;
			plot(t_smpl(:,probe), real(z_p_DM(:,probe)), 'b:', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$p$ [Pa]', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			ylim([-10000,10000]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'NorthEast');
			end
		end%for probe
	end%(plot_outlet_spectra)


	if (plot_SPL)
		g1 = figure();
		set(g1, 'Position', [0 0 1500 250]);
		g2 = figure();
		set(g2, 'Position', [0 0 1500 250]);
		g3 = figure();
		set(g3, 'Position', [0 0 1500 250]);
		for probe = 1:Nprobes
			figure(g1);
			subplot(1,8,probe);
			plot(omega_smpl, p_SPL_MC(:,probe), 'k-', 'LineWidth', 2);
			hold on;
			plot(omega_smpl, p_SPL_DM(:,probe), 'b-', 'LineWidth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('SPL [dB]');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([100, 160]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'NorthEast');
			end

			figure(g2);
			subplot(1,8,probe);
			plot(omega_smpl, s_SPL_MC(:,probe), 'k-', 'LineWidth', 2);
			hold on;
			plot(omega_smpl, s_SPL_DM(:,probe), 'b:', 'LineWidth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('SPL [dB]');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([100, 160]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'NorthEast');
			end

			figure(g3);
			subplot(1,8,probe);
			plot(omega_smpl, z_SPL_MC(:,probe), 'k-', 'lineWidth', 2);
			hold on;
			plot(omega_smpl, z_SPL_DM(:,probe), 'b:', 'LineWIdth', 2);
			xlabel('f [Hz]');
			if (probe == 1)
				ylabel('SPL [dB]');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([-2000, 2000]);
			ylim([100, 160]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'NorthEast');
			end
		end%for probe

		figure(g1);
		subplot(1,8,8);
		plot(omega_smpl, p_SPL_bar_MC, 'k-', 'LineWidth', 2);
		hold on;
		plot(omega_smpl, p_SPL_bar_DM, 'b-', 'LineWidth', 2);
		xlabel('f [Hz]');
		title('Average');
		xlim([-2000, 2000]);
		ylim([100, 160]);
		set(gca, 'FontSize', 16);
		tightfig();

		figure(g2);
		subplot(1,8,8);
		plot(omega_smpl, s_SPL_bar_MC, 'k-', 'LineWidth', 2);
		hold on;
		plot(omega_smpl, s_SPL_bar_DM, 'b-', 'LineWidth', 2);;
		title('Average');
		xlabel('f [Hz]');
		xlim([-2000, 2000]);
		ylim([100, 160]);
		set(gca, 'FontSize', 16);
		tightfig();

		figure(g3);
		subplot(1,8,8);
		plot(omega_smpl, z_SPL_bar_MC, 'k-', 'LineWidth', 2);
		hold on;
		plot(omega_smpl, z_SPL_bar_DM, 'b-', 'LineWidth', 2);
		title('Average');
		xlabel('f [Hz]');
		xlim([-2000, 2000]);
		ylim([100, 160]);
		set(gca, 'FontSize', 16);
		tightfig();
	end%plot_SPL

	if (plot_good_SPL)
		h1 = figure();
		plot(omega_smpl, p_SPL_bar_MC, 'k-', 'lineWidth', 2);
		hold on;
		plot(omega_smpl, p_SPL_bar_DM, 'b:', 'LineWIdth', 2);
		xlabel('f [Hz]');
		ylabel('SPL [dB]');
%		title('Acoustic Forcing');
		xlim([0, 2000]);
		ylim([100, 160]);
		set(gca, 'FontSize', 16);
		h = legend('Compact', 'Non-compact');
		set(h, 'FontSize', 14, 'Location', 'NorthEast');

		h2 = figure();
		plot(omega_smpl, s_SPL_bar_MC, 'k-', 'lineWidth', 2);
		hold on;
		plot(omega_smpl, s_SPL_bar_DM, 'b:', 'LineWIdth', 2);
		xlabel('f [Hz]');
		ylabel('SPL [dB]');
%		title('Entropic Forcing');
		xlim([0, 2000]);
		ylim([100, 160]);
		set(gca, 'FontSize', 16);
		h = legend('Compact', 'Non-compact');
		set(h, 'FontSize', 14, 'Location', 'NorthEast');

		h3 = figure();
		plot(omega_smpl, z_SPL_bar_MC, 'k-', 'lineWidth', 2);
		hold on;
		plot(omega_smpl, z_SPL_bar_DM, 'b:', 'LineWIdth', 2);
		xlabel('f [Hz]');
		ylabel('SPL [dB]');
%		title('Compositional Forcing');
		xlim([0, 2000]);
		ylim([100, 160]);
		set(gca, 'FontSize', 16);
		h = legend('Compact', 'Non-compact');
		set(h, 'FontSize', 14, 'Location', 'NorthEast');

		h4 = figure();
		plot(omega_smpl, t_SPL_bar_MC, 'k-', 'lineWidth', 2);
		hold on;
		plot(omega_smpl, t_SPL_bar_DM, 'b:', 'LineWIdth', 2);
		xlabel('f [Hz]');
		ylabel('SPL [dB]');
%		title('Total Forcing');
		xlim([0, 2000]);
		ylim([100, 160]);
		set(gca, 'FontSize', 16);
		h = legend('Compact', 'Non-compact');
		set(h, 'FontSize', 14, 'Location', 'NorthEast');
	end

	if (plot_inlet_signal)
		f1 = figure();
		set(f1, 'Position', [0 0 1500 250]);
		f2 = figure();
		set(f2, 'Position', [0 0 1500 250]);
		f3 = figure();
		set(f3, 'Position', [0 0 1500 250]);
		for probe = 1:Nprobes

			pi_plus1(:,probe) = fftshift(fft(pi_plus1(:,probe)));
			sigma1(:,probe) = fftshift(fft(sigma1(:,probe)));
			xi1(:,probe) = fftshift(fft(xi1(:,probe)));
			for i = 1:N
				if (abs(He(i)) < 2)
				else
					pi_plus1(i,probe) = 0;
					sigma1(i,probe) = 0;
					xi1(i,probe) = 0;
				end
			end
			pi_plus1(:,probe) = ifft(fftshift(pi_plus1(:,probe)));
			sigma1(:,probe) = ifft(fftshift(sigma1(:,probe)));
			xi1(:,probe) = ifft(fftshift(xi1(:,probe)));

			figure(f1);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), pi_plus1(:,probe), 'k-', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$p$ [Pa]', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end

			figure(f2);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), sigma1(:,probe), 'k-', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$p$ [Pa]', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end
			
			figure(f3);
			subplot(1,7,probe);
			plot(t_smpl(:,probe), xi1(:,probe), 'k-', 'LineWidth', 2);
			xlabel('t [ms]');
			if (probe == 1)
				ylabel('$p$ [Pa]', 'Interpreter', 'LaTeX');
			end
			title(strcat('Sample #', num2str(probe)));
			xlim([min(t_smpl(:,probe)), max(t_smpl(:,probe))]);
			set(gca, 'FontSize', 16);
			if (probe == 1)
				h = legend('Compact', 'Non-compact');
				set(h, 'FontSize', 14, 'Location', 'South');
			end
		end%for probe
	end%(plot_inlet_spectra)

end%function
