#!/usr/bin/python3

import itertools
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.lines import Line2D


import plato_pylib.plato.parse_gau_files as parseGau
import plato_pylib.plato.parse_bas_files as parseBas
import plato_pylib.plato.parse_tbint_files as parseTbint




def gauVsBasPlotSingleFit(gauFile,basFile, field, fieldIdx=0):
	''' Valid field values are nlpp,orbitals, neutatom or density '''
	fittedData = getFieldDatafromGauFile(gauFile, field, fieldIdx)
	outFig = gauVsBasPlotSingleFit_fromGauPolyBas(fittedData, basFile, field, fieldIdx)
	return outFig



def gauVsBasPlotSingleFit_fromGauPolyBas(gauPolyBas, basFile, field, fieldIdx=0):
	origData = getFieldDataFromBasFile(basFile, field, fieldIdx)
	outFig = gauVsBasPlotSingleFit_fromBasTupleAndGauPolyBas(gauPolyBas, origData)	
	return outFig


def gauVsBasPlotSingleFit_fromBasTupleAndGauPolyBas(gauPolyBas, basGridData, multFitByRl=False, lVal=None):

	#Need to generate grid data from the Gaussians expansions
	origData = np.array( basGridData )
	fitVals = np.array( origData )
	fitFunct = gauPolyBas.toGaussFunct()

	for rIdx in range(np.shape(fitVals)[0]):
		fitVals[rIdx,1] = fitFunct( fitVals[rIdx,0] )

	#For pseudopot case we may need to apply r^l factor to term of the fit
	if multFitByRl:
		for rIdx in range(np.shape(fitVals)[0]):
			fitVals[rIdx,1] *= fitVals[rIdx,0]**lVal

	outFig = makePlotFitVsActual(origData,fitVals)

	return outFig

def getFieldDatafromGauFile(gauFile,field,fieldIdx=0):
	parsedGau = parseGau.parseGauFile(gauFile)
	dataSection = parsedGau[field.lower()]
	#orbitals/nl-pp are lists while density/potential are single objs.
	try:
		fieldData = dataSection[fieldIdx]
	except TypeError:
		fieldData = dataSection
	return fieldData


def getFieldDataFromBasFile(basFile,field,fieldIdx=0):
	parsedBas = parseBas.parseBasFile(basFile)
	dataSection = parsedBas[field.lower()] 
	try:
		fieldData = dataSection[fieldIdx]
	except TypeError:
		fieldData = dataSection

	if field.lower() == "orbitals":
		xyData = fieldData.getGridValsPlatoAOrep()
	else:
		xyData = fieldData.gridVals

	return xyData



def makePlotFitVsActual(origXY, fitXY, **kwargs):
	figure = plt.figure()
	ax1 = figure.add_subplot(111)

	ax1.plot( origXY[:,0] , origXY[:,1], label="Original")
	ax1.plot( fitXY[:,0], fitXY[:,1], label="Fit")

	plt.xlabel("Distance")

	return figure


def getFitErrorGaussFit(fitParams:"GauPolyBas obj", origDataGrid:"list of (x,y) tuples for orig data", multFitByRl=False, lVal=None):
	origData = np.array( origDataGrid )
	fitVals = np.array( origData )
	fitFunct = fitParams.toGaussFunct()

	for rIdx in range(np.shape(fitVals)[0]):
		fitVals[rIdx,1] = fitFunct( fitVals[rIdx,0] )

	#For pseudopot case we may need to apply r^l factor to term of the fit
	if multFitByRl:
		for rIdx in range(np.shape(fitVals)[0]):
			fitVals[rIdx,1] *= fitVals[rIdx,0]**lVal

	#Get sum of sqr resids:
	sumSqrRes = sum([ (act-fit)**2 for act,fit in itertools.zip_longest(origData[:,1],fitVals[:,1]) ])

	return sumSqrRes






def plotEosFitVsActual(propDataDict, **kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}
	title = kwargs.get("title", None)

	#Extract data into nice format
	xVals, yVals = list(), list()
	fitX, fitY = list(), list()
	for x,y in propDataDict["data"]:
		xVals.append(x), yVals.append(y)
	for x,y in propDataDict["fitdata"]:
		fitX.append(x), fitY.append(y)

	figure = plt.figure()
	ax1 = figure.add_subplot(111)

	ax1.scatter(xVals,yVals,label="calculated")
	ax1.plot(fitX,fitY,label="fit")
	
	if title is not None:
		plt.title(title)
	plt.xlabel("Volume per atom / $bohr^{3}$")
	plt.ylabel("Energy per atom / eV")
	plt.legend()


	return figure



def createBulkModPlotsFromPropDict(propDicts:"list of dicts",labelList, **kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}

	addLines = kwargs.get("addlines",False)
	title = kwargs.get("title",None)

	#Extract data
	xDataCalcList, yDataCalcList = list(), list()
	for pDict in propDicts:
		currXData, currYData = list(), list()
		for x,y in pDict["data"]:
			currXData.append(x), currYData.append(y)
		xDataCalcList.append(currXData), yDataCalcList.append(currYData)
		
	#Create Figure
	figure = plt.figure()
	ax1 = figure.add_subplot(111)
	plt.rc('text', usetex=True)
	colorOrder = ['darkorange','green','blue', 'red', 'black', 'dimgray','darkturquoise']

	for idx, (xData,yData) in enumerate(zip(xDataCalcList, yDataCalcList)):
		if addLines:
			ax1.plot( xData, yData,  '-o' ,color=colorOrder[idx%len(colorOrder)], label=labelList[idx],)
		else:
			ax1.scatter( xData, yData, color=colorOrder[idx%len(colorOrder)], label=labelList[idx] )
	plt.legend()
	plt.xlabel("Volume per atom / $bohr^{3}$")
	plt.ylabel("Energy per atom / eV")
	
	if title is not None:
		plt.title(title)

	return figure


def createEnergyVsVolCurves_multiStructEachMethod(plotData:"list of lists, first idx is method, second is for diff structures", **kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}

	title = kwargs.get("title",None)
	structLabels = kwargs.get("structlabels", None)
	modelLabels = kwargs.get("modellabels", None)
	xLabel = kwargs.get("xlabel",None)
	yLabel = kwargs.get("ylabel",None)
	xLim = kwargs.get("xlim",None)
	yLim = kwargs.get("ylim",None)

	figure = plt.figure()
	ax1 = figure.add_subplot(111)
	plt.rc('text', usetex=True)


	methodLineStyles = ['-',':']
	structColors = ['darkorange','green','blue', 'red', 'black', 'dimgray','darkturquoise']
	markerStyle = 'x'
	for mIdx, methodData in enumerate(plotData):
		for sIdx,structData in enumerate(methodData):
			ax1.plot(structData[:,0],structData[:,1], lineStyle=methodLineStyles[mIdx], marker=markerStyle, color=structColors[sIdx])

	if title is not None:
		plt.title(title)
	if (modelLabels is not None) and (structLabels is not None):
		legLabels = _getLegNamesEvolCurves(structLabels, modelLabels)
		plt.legend(legLabels)
	if xLabel is not None:
		plt.xlabel(xLabel)
	if yLabel is not None:
		plt.ylabel(yLabel)
	if xLim is not None:
		ax1.set_xlim(xLim[0],xLim[1])
	if yLim is not None:
		ax1.set_ylim(yLim[0],yLim[1])


	return figure
	


def _getLegNamesEvolCurves(structLabels, modelLabels):
	allLabels = list()
	for modLabel in modelLabels:
		for structName in structLabels:
			allLabels.append( "{} , {}".format(modLabel, structName) )

	return allLabels




def plotDos( dosData:"list of nx2 np arrays", **kwargs):
    kwargs = {k.lower():v for k,v in kwargs.items()}

    labels = kwargs.get('labels',None)
    xlim = kwargs.get('xlim',None)
    ylim = kwargs.get('ylim',None)
    title = kwargs.get('title',None)

    figure = plt.figure()
    ax1 = figure.add_subplot(111)
    lineStyles = ['-','--',':']
    
    if labels is None:
        for data in dosData:
            ax1.plot( data[:,0], data[:,1] )
    else:
        for data,currLabel in itertools.zip_longest(dosData,labels):
            ax1.plot( data[:,0], data[:,1], label=currLabel )

        
    #Set the line-styles etc of the data
    dSeries = ax1.lines
    for idx,currLine in enumerate(dSeries):
        currIdx = idx % len(lineStyles)
        currLine.set_linestyle( lineStyles[currIdx] )
    
    plt.legend()
    plt.xlabel("Energy / eV")
    plt.ylabel("Density of States (DoS)")
    if xlim is not None:
        ax1.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax1.set_ylim(ylim[0],ylim[1])
    if title is not None:
        plt.title(title)
    
    return figure






def _plotAllInvSkResults(dataListDict, shellIdxToAngMomA, shellIdxToAngMomB, **kwargs):
    kwargs = {k.lower():v for k,v in kwargs.items()}
    refData = kwargs.get("refData".lower(),None)
    shareX = kwargs.get("sharex","col")
    shareY = kwargs.get("sharey","row")
    xlimAll = kwargs.get("xlim",None)
    ylimAll = kwargs.get("ylim",None)
    xLabel = kwargs.get("xLabel".lower(),None)
    yLabel = kwargs.get("yLabel".lower(),None)
    title = kwargs.get("title",None)
    
    nShellsA, nShellsB = len(dataListDict), len(dataListDict[0])
    angMomToLetter = {0:"s",1:"p",2:"d"}

    nRows, nCols =2, 3
    outFig, allAxes = plt.subplots(nRows, nCols, figsize=(10,7), sharex=shareX, sharey=shareY)
    nPlotted = 0
    currRow, currCol = 0,0
    invSKColors = {"sigma":'b', "pi":'orange', "delta":'green'}
    tbintColors = {"sigma": 'turquoise', "pi":'y', "delta":'limegreen'} 
    markerStyle = 'o'
    #For now just plot all the sigma
    for sIdxA in range(nShellsA):
        for sIdxB in range(nShellsB):
            angMomStr = angMomToLetter[shellIdxToAngMomA[sIdxA]] +  angMomToLetter[shellIdxToAngMomB[sIdxB]]
            legend_elements, legLabels = list(), list()
            if sIdxA <= sIdxB:
                if refData is not None:
                    allAxes[currRow,currCol].plot( refData[sIdxA][sIdxB]["sigma"][:,0], refData[sIdxA][sIdxB]["sigma"][:,1], c=tbintColors["sigma"] )
                allAxes[currRow,currCol].plot( dataListDict[sIdxA][sIdxB]["sigma"][:,0], dataListDict[sIdxA][sIdxB]["sigma"][:,1],
                                              ms=1, marker=markerStyle, linestyle='', c=invSKColors["sigma"])
                legend_elements.append( Line2D([0], [0], linestyle='None',marker=markerStyle, color=invSKColors["sigma"],
                                                markerfacecolor=invSKColors["sigma"], markersize=2) )
                legLabels.append(angMomStr + "$\sigma$")
                
                if dataListDict[sIdxA][sIdxB]["pi"] is not None:
                    currBond = "pi"
                    if refData is not None:
                        allAxes[currRow,currCol].plot( refData[sIdxA][sIdxB][currBond][:,0], refData[sIdxA][sIdxB][currBond][:,1], c=tbintColors[currBond] )
                    allAxes[currRow,currCol].plot( dataListDict[sIdxA][sIdxB][currBond][:,0], dataListDict[sIdxA][sIdxB][currBond][:,1],
                                                    ms=1, marker=markerStyle, linestyle='', c=invSKColors[currBond])
                    legend_elements.append( Line2D([0], [0], linestyle='None',marker=markerStyle, color=invSKColors[currBond],
                                                   markerfacecolor=invSKColors[currBond], markersize=2) )
                    legLabels.append(angMomStr + "$\pi$")

                if dataListDict[sIdxA][sIdxB]["delta"] is not None:
                    currBond = "delta"
                    if refData is not None:
                        allAxes[currRow,currCol].plot( refData[sIdxA][sIdxB][currBond][:,0], refData[sIdxA][sIdxB][currBond][:,1], c=tbintColors[currBond] )
                    allAxes[currRow,currCol].plot( dataListDict[sIdxA][sIdxB][currBond][:,0], dataListDict[sIdxA][sIdxB][currBond][:,1],
                                                    ms=1, marker=markerStyle, linestyle='', c=invSKColors[currBond])
                    legend_elements.append( Line2D([0], [0], linestyle='None',marker=markerStyle, color=invSKColors[currBond],
                                                   markerfacecolor=invSKColors[currBond], markersize=2) ) 
                    legLabels.append(angMomStr + "$\delta$")
                
                allAxes[currRow,currCol].legend(legend_elements,legLabels)
                
                nPlotted += 1
                if currCol + 1 >= nCols:
                    currCol = 0
                    currRow +=1
                else:
                    currCol += 1
                    
                    
    for a in allAxes.flat:
        if xLabel is not None:
            a.set(xlabel=xLabel)
        if yLabel is not None:
            a.set(ylabel=yLabel)
        if xlimAll is not None:
            a.set(xlim=xlimAll)
        if ylimAll is not None:
            a.set(ylim=ylimAll)
    for a in allAxes.flat:
        a.label_outer()
    
    if title is not None:
        allAxes[0,0].set_title(title)
    
    outFig.subplots_adjust(hspace=0.05, wspace=0.05)                    
    return outFig 





if __name__ == '__main__':
	import os
	basePath = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/opt_basis/10el_PP_dorbs/att6_longer_range/create_basis_sets/rc_5pt5/work_folder"
	basPath = os.path.join(basePath,"Mg_denfit.bas")
	gauPath = os.path.join(basePath,"Mg_denfit.gau")
	field = "density"
	currFig = gauVsBasPlotSingleFit(gauPath,basPath,field)
#	makePlotFitVsActual(None,None)

