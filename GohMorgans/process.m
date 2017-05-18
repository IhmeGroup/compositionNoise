function[] = process()
	close all;
	data = load('GohMorgansTransonicProfile.mat');
	x = data.GohMorgansTransonicProfile(:,1);
	M = data.GohMorgansTransonicProfile(:,2);
	plot(x,M);

	AoverAstar = 1./M.*(2/(1.4+1)*(1+0.4/2.*M.*M)).^(2.4/(2*.4));

	plot(x, AoverAstar);
end
