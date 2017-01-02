function[phase] = unwrap(phase)
	dphase = diff(phase);
	for i = 1:length(dphase)
		if (abs(dphase(i)) > pi)
			phase(i+1:end) = phase(i+1:end) - 2*pi;
			dphase = diff(phase);
		end
	end
end
