

import itertools as it

import gen_basis_helpers.shared.surfaces as surfHelp
import gen_basis_helpers.adsorption.water_adsorbate as waterAdsHelp


def getStructWithNBilayersAdded_standard(inpCell, nLayers, topLayerWaterAdsDetector,
										 bilayerSpacing, layerTol, absVacLen=None, topAdsTemplate=None):
	""" Takes a water bilayer structure returns a structure with n-added bilayers
	
	Args:
		inpCell: (plato_pylib UnitCell object) The starting (water bilayer) structure
		nLayers: (int) Number of layers to add (must be >=1)
		topLayerWaterAdsDetector: (adsorption.detect_water_layer.DetectOuterAdsorbedWaterLayer object) 
		bilayerSpacing: (float) The closest O-O gap between two bilayers
		layerTol: (float) Two O will be considered in the same "layer" if they are closer than this
		absVacLen: (Optional, float) If None then output will have the same length as input cell. If absVacLen is set then this is the vacuum length in the cell
		topAdsTemplate: (Optional, WaterAdsorbateStandard) If None then we use the "top" molecules of the current bilayer as a template for building the others. Idea is to set this to Hup or Hdown to build diff bilayers

	Returns
		outCell: (plato_pylib UnitCell object) With the additional bilayers added
 
	"""

	#0) Simple edge-case handling
	if (nLayers < 1):
		raise ValueError("nLayers must be 1 or greater; {} is not allowed".format(nLayers))

		
	#4) Create the output cell; 
	estLayerHeight = bilayerSpacing + topLayerWaterAdsDetector.waterIdxDetector.maxHeightLayer
	estNewLayersHeight = estLayerHeight*nLayers
	outCell = surfHelp.GenericSurface(inpCell, 1, lenVac=estNewLayersHeight*2).unitCell

		
	#1) Get the current top-adsorbed layer
	startLayerAdsObjs = topLayerWaterAdsDetector.getAdsorbateObjsFromInpGeom(outCell)

	#2) Create the first extra layer (using topLayerTemplate)
	args = [startLayerAdsObjs, bilayerSpacing]
	kwargDict = {"layerTolerance":layerTol, "surfNormal":[0,0,1], "topLayerTemplate":topAdsTemplate}
	extraLayers = [waterAdsHelp.getAdsorbateObjsForNextWaterBilayerBasic(*args,**kwargDict)]

	#3) Create any extra layers
	kwargDict.pop("topLayerTemplate")
	for idx in range(1,nLayers):
		args[0] = extraLayers[-1]
		currLayer = waterAdsHelp.getAdsorbateObjsForNextWaterBilayerBasic(*args,**kwargDict)
		extraLayers.append(currLayer)


	# #5) Add co-ords to output cell AND set the absolute vacuum region
	newCoords = list()
	for adsObj in it.chain(*extraLayers):
		newCoords.extend(adsObj.geom)
	outCell.cartCoords += newCoords

	if absVacLen is not None:
		outCell = surfHelp.GenericSurface(outCell, 1, lenAbsoluteVacuum=absVacLen).unitCell
	else:
		outCell = surfHelp.GenericSurface(outCell, 1, lenVac=0).unitCell

	return outCell




