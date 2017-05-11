function[p_out, t_out] = cleanProbes(t_in, p_in, probe_cleaning)

	if (probe_cleaning)
		figure();
		plot(t_in, p_in, 'LineWidth', 2);
		xlabel('t [ms]');
		ylabel('\phi [??]');
		title('Raw Signal');
		set(gca, 'FontSize', 16);
	end

	p_in = p_in(1:end-4);
	N = length(p_in);
	M = N/4;
	t_1 = t_in(0*M+1:1*M);
	p_1 = p_in(0*M+1:1*M);
	t_2 = t_in(1*M+1:2*M);
	p_2 = p_in(1*M+1:2*M);
	t_3 = t_in(2*M+1:3*M);
	p_3 = p_in(2*M+1:3*M);
	t_4 = t_in(3*M+1:4*M);
	p_4 = p_in(3*M+1:4*M);
	t_5 = t_in(0.5*M+1:1.5*M);
	p_5 = p_in(0.5*M+1:1.5*M);
	t_6 = t_in(1.5*M+1:2.5*M);
	p_6 = p_in(1.5*M+1:2.5*M);
	t_7 = t_in(2.5*M+1:3.5*M);
	p_7 = p_in(2.5*M+1:3.5*M);

	p_1 = p_1 - mean(p_1);
	p_2 = p_2 - mean(p_2);
	p_3 = p_3 - mean(p_3);
	p_4 = p_4 - mean(p_4);
	p_5 = p_5 - mean(p_5);
	p_6 = p_6 - mean(p_6);
	p_7 = p_7 - mean(p_7);

	h = hanning(M);
	area = trapz(h)./M;
	h = h./area;

	p_1 = p_1.*h;
	p_2 = p_2.*h;
	p_3 = p_3.*h;
	p_4 = p_4.*h;
	p_5 = p_5.*h;
	p_6 = p_6.*h;
	p_7 = p_7.*h;

	p_out = [p_1 p_2 p_3 p_4 p_5 p_6 p_7];
	t_out = [t_1 t_2 t_3 t_4 t_5 t_6 t_7];
	if (probe_cleaning)
		h = figure();
		set(h, 'Position', [0 0 1000, 400]);
		subplot(2,4,1);
		plot(t_1, p_1, 'LineWidth', 2);
		xlabel('t [ms]');
		ylabel('\phi [??]');
		title('Sample 1');
		set(gca, 'FontSize', 14);
		xlim([min(t_1), max(t_1)])

		subplot(2,4,2);
		plot(t_2, p_2, 'LineWidth', 2);
		xlabel('t [ms]');
		ylabel('\phi [??]');
		title('Sample 2');
		set(gca, 'FontSize', 14);
		xlim([min(t_2), max(t_2)])

		subplot(2,4,3);
		plot(t_3, p_3, 'LineWidth', 2);
		xlabel('t [ms]');
		ylabel('\phi [??]');
		title('Sample 3');
		set(gca, 'FontSize', 14);
		xlim([min(t_3), max(t_3)])

		subplot(2,4,4);
		plot(t_4, p_4, 'LineWidth', 2);
		xlabel('t [ms]');
		ylabel('\phi [??]');
		title('Sample 4');
		set(gca, 'FontSize', 14);
		xlim([min(t_4), max(t_4)])

		subplot(2,4,5);
		plot(t_5, p_5, 'LineWidth', 2);
		xlabel('t [ms]');
		ylabel('\phi [??]');
		title('Sample 5');
		set(gca, 'FontSize', 14);
		xlim([min(t_5), max(t_5)])

		subplot(2,4,6);
		plot(t_6, p_6, 'LineWidth', 2);
		xlabel('t [ms]');
		ylabel('\phi [??]');
		title('Sample 6');
		set(gca, 'FontSize', 14);
		xlim([min(t_6), max(t_6)])

		subplot(2,4,7);
		plot(t_7, p_7, 'LineWidth', 2);
		xlabel('t [ms]');
		ylabel('\phi [??]');
		title('Sample 7');
		set(gca, 'FontSize', 14);
		xlim([min(t_7), max(t_7)])

	end

end
