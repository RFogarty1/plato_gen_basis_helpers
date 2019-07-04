
""" Purpose of this code is to help me run calculations to get the bulk-modulus using CP2K """

import os
import pathlib

import bulk_mod_fitting as fitBMod
import cp2k_file_helpers as cp2kHelp
import cp2k_calc_objs as cp2kCalcObjs


#Specific verison of the CalcObj which we need for this
CP2KCalcObj = cp2kCalcObjs.CP2KCalcObj
cp2kCalcObjs.addInpPathDescriptorToCP2KCalcObjCLASS(CP2KCalcObj)



def getBulkModDictsFromCalcObjDict(calcObjsDict):
    allBModDicts = dict()
    for key, currObjs in calcObjsDict.items():
        currOutFiles = [x.outFilePath for x in currObjs]
        currBModDict = fitBMod.getBulkModFromOutFilesAseWrapper(currOutFiles)
        allBModDicts[key] = currBModDict
    return allBModDicts



def createObjsDictFromOptions(startFolder, structsDict, elementBasisInfo, basisFile, potFile, kptsDict,
                             absGridValsDict, relGridValsDict):
    outObjsDict = dict()
    for key,currStruct in structsDict.items():
        currFolder = os.path.abspath(os.path.join(startFolder,key))
        currKPts = kptsDict[key]
        absGridVals, relGridVals = absGridValsDict[key], relGridValsDict[key]
        currObj = _createBulkModFilesOneStruct(currFolder, currStruct, elementBasisInfo, basisFile, potFile, currKPts,
                                              absGridVals, relGridVals)
        outObjsDict[key] = currObj
    return outObjsDict



def _createBulkModFilesOneStruct(folderPath:str, uCells:list, elementBasisInfo, basisFile, potFile, kPts:list, absCutoff:float, relCutoff:float):
    pathlib.Path(folderPath).mkdir(parents=True,exist_ok=True)
    
    allFileObjs = list()
    for currStruct in uCells:
        currObj = _createCP2KSinglePointEnergyCalcObj(folderPath, currStruct, elementBasisInfo, basisFile, potFile, kPts, absCutoff, relCutoff)
        allFileObjs.append(currObj)
    
    for x in allFileObjs:
        x.writeFile()
    
    return allFileObjs

def _createCP2KSinglePointEnergyCalcObj(folderPath, uCell, elementBasisInfo, basisFile, potFile, kPts, absCutoff, relCutoff):
    modDictAll = {"kpts":kPts,
                 "outdir":folderPath,
                 "basisFile":basisFile,
                 "potFile":potFile,
                 "gridCutRel":relCutoff,
                 "gridCutAbs":absCutoff}
    
    nameFmt = "bmod_calcs_a_{}_relgrid{}_absgrid{}"
    aValue = "{:.4f}".format( uCell.lattParams["a"] ).replace(".","pt")
    outFileName = nameFmt.format(aValue, relCutoff, absCutoff)  
    outCp2kObj = cp2kHelp.createCp2kSinglePointObjFromUCellGeom(uCell, elementBasisInfo, **modDictAll)
    
    basePath = os.path.join(folderPath,outFileName)
    outObj = CP2KCalcObj(outCp2kObj,basePath=basePath)


    return outObj



def getRunCommsFromDictOfCalcObjLists(calcObjsDict):
    allRunComms = list()
    for currObjs in calcObjsDict.values(): #Dont care about the order
        currComms = [x.runComm for x in currObjs]
        allRunComms.extend(currComms)
    return allRunComms



def printAllEosFitResults(eqnStateDicts, structOrder=None):
    allE0Vals = [x["e0"] for x in eqnStateDicts.values()]
    allDeltaE0Vals = { key:val["e0"]-min(allE0Vals) for key,val in eqnStateDicts.items() }

    if structOrder is None:
        structOrder = eqnStateDicts.keys() #Print in a random order
    
    print("Structure, V0, B0, e0")
    
    for label in structOrder:
        fitRes, deltaE0 = eqnStateDicts[label], allDeltaE0Vals[label]
        print("{}, {:.3f}, {:.3f}, {:.3f}".format(label, fitRes["v0"], fitRes["bulkMod"], deltaE0))
    

