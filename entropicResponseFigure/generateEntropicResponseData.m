function[] = generateEntropicResponseData()
	addpath('../HelmholtzSweep');
	[OMEGA, TRANS, PHASE] = generateHelmholtzSweepData(3);
	save('entropicResponseData.mat', 'OMEGA', 'TRANS', 'PHASE');
end
