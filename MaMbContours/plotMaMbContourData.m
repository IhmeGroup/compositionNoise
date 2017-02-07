function[] = plotMaMbContourData()
	close all;

	h = dir;
	N = length(h);

	for file = 1:N
		filename = h(file).name;
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
end
