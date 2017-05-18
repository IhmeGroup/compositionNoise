function[M] = MFromEtaDMSC(eta, SPLINES);
	gamma = returnAmbientState();
	gp1 = gamma + 1;
	gm1 = gamma - 1;
	gp1o2 = gp1/2;
	gm1o2 = gm1/2;
	if ((~exist('DMSC1')) || (~exist('DMSC2')))
		data = load('../GohMorgans/GohMorgansShockProfile.mat');
		ETA = data.GohMorgansShockProfile(:,1);
		MACH = data.GohMorgansShockProfile(:,2);
		ETA1 = smooth(smooth(ETA(1:49)));
		MACH1 = smooth(smooth(MACH(1:49)));
		ETA2 = smooth(smooth(ETA(50:end)));
		MACH2 = smooth(smooth(MACH(50:end)));
		global DMSC1;
		global DMSC2;
		DMSC1 = spline(ETA1, MACH1);
		DMSC2 = spline(ETA2, MACH2);
	end
	if (eta <= 0.5008)
		M = ppval(DMSC1, eta);
	else
		M = ppval(DMSC2,eta);
	end
end
