function[] = generateAcousticResponseData()
	[OMEGA, TRANS, PHASE] = generateJeonglaeSweepData(1);
	save('acousticResponseData.mat', 'OMEGA', 'TRANS', 'PHASE');
end
