function[] = plotMaMbContourData()
	close all;

	fs = 18;


	filenames = {	'MaMbContoursData.0.mat';
					'MaMbContoursData.0.5.mat';
					'MaMbContoursData.1.mat';
					'MaMbContoursData.1.5.mat';
					'MaMbContoursData.2.mat';};

	f = figure();
	set(f, 'Position', [0 0 1200 600]);
	g = figure();
	set(g, 'Position', [0 0 1200 600]);
	h = figure();
	set(h, 'Position', [0 0 1200 600]);
	for file = 1:5
		filename = char(filenames(file))
		load(filename);

		[N,M] = size(MA);
		for i = 1:N
			for j = M:-1:1
				if (MA(i,j) == MB(i,j))
					WPSabove = (WPS(i,j+2) - WPS(i,j+1))/(MB(i,j+2) - MB(i,j+1))*(MB(i,j+1) - MB(i,j)) + WPS(i,j+1);
					WPZabove = (WPZ(i,j+2) - WPZ(i,j+1))/(MB(i,j+2) - MB(i,j+1))*(MB(i,j+1) - MB(i,j)) + WPZ(i,j+1);
					if (i > 2)
						WPSleft = (WPS(i-2,j) - WPS(i-1,j))/(MA(i-2,j) - MA(i-1,j))*(MA(i-1,j) - MA(i,j)) + WPS(i-1,j);
						WPZleft = (WPZ(i-2,j) - WPZ(i-1,j))/(MA(i-2,j) - MA(i-1,j))*(MA(i-1,j) - MA(i,j)) + WPZ(i-1,j);
					else
						WPSleft = WPSabove;
						WPZleft = WPZabove;
					end
					WPS(i,j) = 0.5*(WPSabove + WPSleft);
					WPZ(i,j) = 0.5*(WPZabove + WPZleft);
				end
			end
		end




		for run = 1:3
			if (run == 1)
				figure(f);
			elseif (run == 2)
				figure(g);
			elseif (run == 3)
				figure(h);
			end
			if (file == 1)
				subplot('Position', [0.05 0.55 .8/3 0.4])
			elseif (file == 2)
				subplot('Position', [0.05*2+.8/3 0.55 .8/3 0.4])
			elseif (file == 3)
				subplot('Position', [0.05*3+2*.8/3 0.55 .8/3 0.4])
			elseif (file == 4)
				subplot('Position', [.2083 0.05 .8/3 .4])
			elseif (file == 5)
				subplot('Position', [.2083+.8/3+0.05 0.05 .8/3 .4])
			end
			if (run == 1)
				surface(MA, MB, abs(WPZ) + 1E-3, 'EdgeColor', 'none');
			elseif (run == 2)
				surface(MA, MB, abs(WPS) + 1E-3, 'EdgeColor', 'none');
			elseif (run == 3)
				surface(MA, MB, abs(WPZ + 1E-3)./abs(WPS + 1E-3), 'EdgeColor', 'none');
			end
			shading interp;
			hold on;
			r = caxis();
			a = fill3([0 1 1], [0 0 1], [100 100 100], 'w');
			set(a, 'EdgeColor','none');
			caxis(r);

			xlabel('$M_a$', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
			ylabel('$M_b$', 'Interpreter', 'LaTeX', 'FontSize', fs, 'FontName', 'Times');
			set(gca, 'FontSize', fs, 'FontName', 'Times');
			if (file == 1) 
				title('$(a)$', 'FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
			elseif (file == 2)
				title('$(b)$', 'FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
			elseif (file == 3)
				title('$(c)$', 'FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
			elseif (file == 4)
				title('$(d)$', 'FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
			elseif (file == 5)
				title('$(e)$', 'FontSize', fs, 'FontName', 'Times', 'Interpreter', 'LaTeX');
			end
			colorbar();
			colormap('jet');
		end%for run
	end%for file

	figure(f);
	print -depsc CompositionContours.eps

	figure(g);
	print -depsc EntropyContours.eps

	figure(h);
	print -depsc RatioContours.eps

end%function

function[] = shitty()

	for run = 1:3				
	for file = 1:1
		filename = char(filenames(file))
		load(filename);
			if (run == 1) 
				figure(f);
			elseif (run == 2)
				figure(g);
			elseif (run == 3)
				figure(h);
			end
			if (i == 1)
				subplot('Position', [0.05 0.55 .8/3 0.4])
			elseif (i == 2)
				subplot('Position', [0.05*2+.8/3 0.55 .8/3 0.4])
			elseif (i == 3)
				subplot('Position', [0.05*3+2*.8/3 0.55 .8/3 0.4])
			elseif (i == 4)
				subplot('Position', [.2083 0.05 .8/3 .4])
			elseif (i == 5)
				subplot('Position', [.2083+.8/3+0.05 0.05 .8/3 .4])
			end
			if (run == 1)
				surface(MA, MB, abs(WPZ) + 1E-3, 'EdgeColor', 'none');
			elseif (run == 2)
				surface(MA, MB, abs(WPS) + 1E-3, 'EdgeColor', 'none');
			elseif (run == 3)
				surface(MA, MB, abs(WPZ + 1E-3)/(WPS + 1E-3), 'EdgeColor', 'none');
			end
			xlabel('$M_a$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_b$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
			colorbar();
			colormap('jet');
		end%for run
	end%for file
end%function[]

function[] = bullshit()




		if strncmpi(filename, 'MaMbContoursData.', 15)
			disp('We got ourselves a winner!')
			load(filename);

			[a,b] = size(WPS);
			for c = 1:a
				for d = 1:b
					if (MA(c,d) > MB(c,d)) 
						WPS(c,d) = Inf; 
						WPZ(c,d) = Inf;
					end
				end
			end

%			Composition Noise Contour Field
			figure();
			surface(MA, MB, abs(WPZ) + 1E-3, 'EdgeColor', 'none');
			shading interp;
			xlabel('$M_a$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_b$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
			colorbar();
			colormap('jet');
			title(strcat('Composition:', filename));

%			Entropy Noise Contour Field
			figure();
			surface(MA, MB, abs(WPS + 1E-3), 'EdgeColor', 'none');
			shading interp;
			xlabel('$M_a$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_b$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
			colorbar();
			colormap('jet');
			title(strcat('Entropy:', filename));
			
%			Ratio of Composition Noise to Entropy Noise Contour Field
			figure();
			surface(MA, MB, abs(WPZ + 1E-3)./abs(WPS + 1E-3), 'EdgeColor', 'none');
			shading interp;
			colorbar();
			xlabel('$M_a$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
			ylabel('$M_b$', 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
			colormap('jet');
			set(gca, 'FontSize', 14, 'FontName', 'Times');
			title(strcat('Ratio: ', filename));
		else
			disp('Sorry, doesnt match bro')
		end%file is a data file

	end%for file
