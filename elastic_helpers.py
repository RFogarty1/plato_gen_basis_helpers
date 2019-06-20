#!/usr/bin/python3

''' Purpose is to provide functions to help calc. elastic constants. These are meant to be higher lvl functs than in plato_pylib '''

import copy
import itertools as it
import os
import pathlib
import evol_test_vs_ref_helpers as evolHelpers
import dft_file_helpers as dftHelp
import calc_methods as calcMethods

import plato_pylib.utils.elastic_consts as elastic
import plato_pylib.plato.mod_plato_inp_files as platoInp
import plato_pylib.parseOther.parse_castep_files as parseCastep
import plato_pylib.plato.parse_plato_out_files as parsePlatoOut


import numpy as np
import matplotlib.pyplot as plt


#def createPlatoElasCalcsFileObjs_fromMethodStr(targFolder:"str, folder path", startStruct:"UnitCell class obj", strainParams:list, crystType, dataSet, extraKwargs:dict)
#	allStrainedStructs = elastic.getStrainedStructsForElasticConsts(startStruct, strainParams, crystType=crystType)
#	allOutObjs = list()
#
#	for idx, structSet in enumerate(allStrainedStructs):
#		currObjs = _createPlatoFilesOneStrainPattern(targFolder, structSet, idx, strainparams, fileStrDict)
#
#
#
#def _createPlatoFilesOneStrainPattern_fromMethodStr(targFolder:"str, folder path", strainedStructs:"list of ucell objs", pattIdx:"int, idx for naming files",
#                                                    strainCoeffs:list, methodStr:str, dataSet, kpts, integGrid=None):
#




def createPlatoElasCalcsFileObjs(targFolder:"str, folder path", startStruct:"UnitCell class obj", strainParams:list, inpOptDict:dict, crystType, dftFiles=False):

	allStrainedStructs = elastic.getStrainedStructsForElasticConsts(startStruct, strainParams, crystType=crystType)

	if dftFiles:
                fileDict = dftHelp.loadDefaultDftOptDict()
	else:
                fileDict = evolHelpers.loadDefaultTb2OptDict()

	fileDict.update(inpOptDict) #User options always overwrite the defaults when given

	if dftFiles:
		fileStrDict = dftHelp.getPlatoStrDictFromOptDict_dft(fileDict)
	else:
		fileStrDict = evolHelpers.getPlatoStrDictFromOptDict_tb1OrTb2(fileDict)

	allOutObjs = list()

	for idx,structSet in enumerate(allStrainedStructs):
		currObjs = _createPlatoFilesOneStrainPattern(targFolder,structSet, idx, strainParams, fileStrDict)
		allOutObjs.append(currObjs)

	return allOutObjs

def _createPlatoFilesOneStrainPattern(targFolder:"str, folder path", strainedStructs:"list of ucell objs", pattIdx:"int, idx for naming files",
                                      strainCoeffs:list, optStrDict:"plato option dict, with all except geom fields defined"):

	#1 create folder + sort out all the file paths
	baseFileFmt = "calc_elastic_{}_{}.in"
	pathlib.Path(targFolder).mkdir(parents=True, exist_ok=True)
	basePath = os.path.abspath( os.path.join(targFolder,baseFileFmt) )

	#Create the files
	outObjs = list()
	for currStrain,currStruct in it.zip_longest(strainCoeffs,strainedStructs):
		currDict = copy.deepcopy(optStrDict) #V.likely overkill, since i only mod geom-fields which should be absent anyway
		geomSection = platoInp.getPlatoGeomDictFromUnitCell(currStruct)
		currDict.update(geomSection)
		strainStr = "{:.4f}".format(currStrain).replace(".","pt").replace("-","m")
		currPath = basePath.format(pattIdx,strainStr)
		platoInp.writePlatoOutFileFromDict(currPath, currDict)
		outObjs.append( ElasticCalcObj(currStrain, pattIdx, currPath) )

	return outObjs

def getStrainCoeffs(minVal:"float, inclusive", maxVal:float, stepSize:float):
    outList = list()
    nVals = (maxVal - minVal) / stepSize
    x = minVal
    while x < (maxVal + stepSize):
        outList.append(x)
        x+=stepSize
    return outList



def getStrainStressFromCalcObjs(objList, outExt:"e.g. .castep"):
	evToJoule = 1.60218e-19
	bohrToMetre = 5.29177e-11
	pascalToGpa = 1/(1e9)
	rydToJoule = 2.1798741E-18
	allPaths = [x.inpPath+outExt for x in objList]


	if outExt == ".castep":
		parser = parseCastep.parseCastepOutfile
		elasticConvFactor = (evToJoule / (bohrToMetre**3))*pascalToGpa
	elif outExt == ".out":
		parser = parsePlatoOut.parsePlatoOutFile
		elasticConvFactor = (rydToJoule / (bohrToMetre**3))*pascalToGpa
	else:
		raise ValueError("{} is an invalid file ext".format(outExt))

	energies, deltaE, volumes = list(), list(), list()
	for x in objList:
		currPath = x.inpPath + outExt
		currParsed = parser(currPath)
		nAtoms = currParsed["numbAtoms"]
		if outExt == ".castep":
			currParsed["unitCell"].convAngToBohr()
		energies.append( currParsed["energies"].electronicTotalE )
		volumes.append( currParsed["unitCell"].volume )

	deltaE = [x-min(energies) for x in energies]
	strains = [x.strainVal for x in objList]
	stresses = [(e/vol)*elasticConvFactor for e,vol in it.zip_longest(deltaE,volumes)]
	strainStresses = [(x,y) for x,y in it.zip_longest(strains,stresses)]
	return strainStresses


def plotElasticSecondDerivFit(fitRes, **kwargs):
	#Parse kwargs
	kwargs = {k.lower():v for k,v in kwargs.items()}
	xLabel = kwargs.get("xlabel", None)
	yLabel = kwargs.get("ylabel", None)
	title = kwargs.get("title", None)
	propDict = {"xlabel": (xLabel, plt.xlabel),
	            "ylabel": (yLabel, plt.ylabel),
	            "title" : (title, plt.title)}


	polyFunct = lambda x: ( fitRes.params[0]*(x**2) ) + ( fitRes.params[1]*x ) + fitRes.params[2]
	data = np.zeros( ( len(fitRes.strainVals),2) )
	
	for idx in range(len(fitRes.strainVals)):
		data[idx,0] = fitRes.strainVals[idx] 
		data[idx,1] = fitRes.stressVals[idx]
		
	#Get the fit data
	fitXVals = np.linspace(data[0,0], data[-1,0])
	fitYVals = [polyFunct(x) for x in fitXVals]
	
	
	#Make the figure
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	
	ax1.scatter(data[:,0],data[:,1], label='calculated')
	ax1.plot(fitXVals,fitYVals, label='fitted')


	#Add stuff
	for key in propDict.keys():
		if propDict[key][0] is not None:
			propDict[key][1](propDict[key][0])


	return fig


class ElasticCalcObj:
    def __init__(self, strainVal, pattIdx,inpPath):
        self.strainVal = strainVal
        self.inpPath = os.path.splitext(inpPath)[0]
        self.pattIdx = pattIdx

