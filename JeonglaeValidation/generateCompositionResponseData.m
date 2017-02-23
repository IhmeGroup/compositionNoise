function[] = generateCompositionResponseData()
	[OMEGA, TRANS, PHASE] = generateJeonglaeSweepData(4);
	save('compositionResponseData.mat', 'OMEGA', 'TRANS', 'PHASE');
end
