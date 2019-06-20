
import evol_test_vs_ref_helpers as evoltest
import plato_pylib.plato.mod_plato_inp_files as platoInp

#DFT_ONLY_OPT_KEYS


def writePlatoDftFileFromDict(outPath,inpDict):
	optDict = loadDefaultDftOptDict()
	optDict = {k.lower():v for k,v in optDict.items()}

	#Update any options requested by usr
	dCopy = {k.lower():v for k,v in inpDict.items()}
	optDict.update(dCopy)

	outDict = getPlatoStrDictFromOptDict_dft(optDict)
	outDict["integralmesh"] = str(0)
	outDict["scfflag"] = str(inpDict["scfflag"])

	platoInp.writePlatoOutFileFromDict(outPath, outDict)


def loadDefaultDftOptDict():
	return platoInp.getDefOptDict("dft")
#	outDict = {k.lower():v for k,v in evoltest.loadDefaultTb2OptDict().items()} #A LOT of overlap between the input files
#
#
#	#DFT specific flags
#	outDict["integralmesh"] = 0
#	outDict["model"] = "lcao" #will correspond to 2
#	outDict["optimizemesh"] = 1
#	outDict["fftGridSpacing"] = 0.7
#	outDict["ACGGridSpacing"] = [30,29,3]
#	outDict["paralleliseFFT"] = 0
#	outDict["SpinMixLevels"] = 5
#	outDict["diagonalisationMethod"] = 1
#	outDict["DiagonalisationWorkSize"] = 500
#	outDict["MaxBond"] = 6.0
#	outDict["DensityRWFlag"] = 0
#	outDict["writeAtomDensityFlag"] = 0
#	outDict["WavefunctionFlag"] = 0
#	outDict["HyperfineFlag"] = 0
#
#	#Parameters that should seldom be chnaged
#	outDict["ACGEwaldParameters"] = [0.4,10.0]
#	outDict["ACGMeshReduction"] = [6, 5.0]
#	outDict["ACGMinPartitionWt"] = 0.0
#	outDict["ACGAngularMeshType"] = 1
#	outDict["DiameterNLV"] = 0.0
#	outDict["IntegralPartitionFlag"] = 0
#	outDict["DensityFit"] = [1,0]
#	outDict["NumericVna"] = 0
#	outDict["SpinMixScheme"] = 1
#	outDict["MixThresHold"] = 100
#	outDict["MixMetric"] = [0, 20.0]
#	outDict["SplitDensityFlag"] = 0
#
#	outDict = {k.lower():v for k,v in outDict.items()}
#
#	return outDict
#

def getPlatoStrDictFromOptDict_dft(inpDict):
	return platoInp.getStrDictFromOptDict(inpDict,"dft")
#	optDict = evoltest.getPlatoStrDictFromOptDict_tb1OrTb2(inpDict)
#
#	modelToFlag = {"lcao":"2"}
#
#	if "XCFunctional".lower() in inpDict.keys():
#		optDict["XCFlag".lower()] = optDict["XCFunctional".lower()]
#
#	if "model" in inpDict.keys():
#		optDict["model"] = modelToFlag[ optDict["model"] ]
#
#	if "ACGGridSpacing".lower() in inpDict.keys():
#		optDict["ACGGridSpacing".lower()] = " ".join([str(x) for x in optDict["ACGGridSpacing".lower()] ])
#
#	if "ACGEwaldParameters".lower() in inpDict.keys():
#		optDict["ACGEwaldParameters".lower()] = " ".join([str(x) for x in optDict["ACGEwaldParameters".lower()] ])
#
#	if "ACGMeshReduction".lower() in inpDict.keys():
#		optDict["ACGMeshReduction".lower()] = " ".join([str(x)for x in optDict["ACGMeshReduction".lower()] ])
#
#	if "DensityFit".lower() in inpDict.keys():
#		optDict["DensityFit".lower()] = " ".join([str(x) for x in optDict["DensityFit".lower() ]])
#
#	if "MixMetric".lower() in inpDict.keys():
#		optDict["MixMetric".lower()] = " ".join([str(x) for x in optDict["MixMetric".lower() ]])
#
#	#Convert any remaining non-strings
#	for k in optDict.keys():
##		if k in evoltest.DFT_ONLY_OPT_KEYS:
#		optDict[k] = str(optDict[k])
#
#
#	optDict["integralmesh"] = str(inpDict["integralmesh"])
#	#"XCFlag" keyword used for the functional in dft
#	#IntegralMesh is same keyword but different fields
#	return optDict
#

