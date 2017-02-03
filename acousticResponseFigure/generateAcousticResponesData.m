function[] = generateAcousticResponseData()
	addpath('../HelmholtzSweep');
	[OMEGA, TRANS, PHASE] = generateHelmholtzSweepData(1);
	save('acousticResponseData.mat', 'OMEGA', 'TRANS', 'PHASE');
end
