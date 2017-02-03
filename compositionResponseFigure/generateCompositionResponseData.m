function[] = generateCompositionResponseData()
	addpath('../HelmholtzSweep');
	[OMEGA, TRANS, PHASE] = generateHelmholtzSweepData(4);
	save('compositionResponseData.mat', 'OMEGA', 'TRANS', 'PHASE');
end
