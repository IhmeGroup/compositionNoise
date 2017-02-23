function[] = generateEntropicResponseData()
	[OMEGA, TRANS, PHASE] = generateJeonglaeSweepData(3);
	save('entropicResponseData.mat', 'OMEGA', 'TRANS', 'PHASE');
end
