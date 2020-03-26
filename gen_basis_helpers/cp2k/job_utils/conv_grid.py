
""" Very specific functions for running CP2K grid-convergence calculations """

import copy
import os
import itertools as it

import gen_basis_helpers.cp2k.cp2k_calc_objs as baseCalcObjCP2K
import gen_basis_helpers.cp2k.basis_register as basRegister
import gen_basis_helpers.cp2k.method_register as methRegister
import gen_basis_helpers.cp2k.cp2k_file_helpers as pyCP2KHelp
import gen_basis_helpers.shared.label_objs as labelHelp
import gen_basis_helpers.shared.calc_runners as calcRunners

import gen_basis_helpers.workflows.convergers as wflowConv


baseCalcObjCP2K.addInpPathDescriptorToCP2KCalcObjCLASS( baseCalcObjCP2K.CP2KCalcObj )
baseCalcObjCP2K.addRelGridDescriptorToCP2KCalcObjCLASS( baseCalcObjCP2K.CP2KCalcObj )
baseCalcObjCP2K.addAbsGridDescriptorToCP2KCalcObjCLASS( baseCalcObjCP2K.CP2KCalcObj )
baseCalcObjCP2K.addMaxScfDescriptorToCP2KCalcObjCLASS ( baseCalcObjCP2K.CP2KCalcObj )


def createGridValsCP2K(varyVals, fixedVal, convType):
	allObjs = [GridValsCP2K(None,None, convType) for x in range(len(varyVals))]
	if convType == "abs":
		varyAttr, fixedAttr = "absCut", "relCut"
	elif convType == "rel":
		varyAttr, fixedAttr = "relCut", "absCut"
	else:
		raise ValueError("gridConvType must be either abs or rel; not {}".format(convType))
	
	for varyVal,currObj in zip(varyVals,allObjs):
		setattr(currObj, varyAttr , varyVal)
		setattr(currObj, fixedAttr, fixedVal)
	
	return allObjs


class GridValsCP2K():
	""" Class used to represent grid cutoff energies in CP2K

	Attributes:
		absCut: Absolute cutoff energy (just called "cutoff" in CP2K)
		relCut: Relative cutoff energy
		varyValCut: The value of the cutoff energy which is being varied in the current grid-based convergence. For example
		            it will return relCut if your converging the relative cutoff
	"""
	def __init__(self, absCut, relCut, varyType):
		""" Initializer
		
		Args:
			absCut: (float) Absolute cutoff energy
			relCut: (float) Relative cutoff energy
			varyType: (str), Type of cutoff your converging. Either "abs" or "rel"
		"""
		self.absCut = absCut
		self.relCut = relCut
		self.varyType = varyType
		
	@property
	def varyValCut(self):
		if self.varyType == "abs":
			return self.absCut
		elif self.varyType == "rel":
			return self.relCut
		else:
			raise ValueError("{} is an invalid varyType".format(varyType))

	def toStr(self):
		return "abs{}_rel{}".format(self.absCut, self.relCut).replace(".","pt")



#TODO: Use context manager to temporarily set attribute of object when calling createStandardInput
#TODO: Probably enforce some arguments and/or write some test code to check
class CreatorForStandardInputForGridConvergenceTemplate():
	"""Class used to create standard input object for running CP2K grid-convergence for multiple basis sets.
	
	createStandardInput is the function that should be called to create the input object

	"""


	def __init__(self, geom=None, kPts=None, maxScf=None, cp2kMethodStr=None,
				basisKeyDicts=None, basisAliases=None, workFolder=None, workFlowToOutputObjMap=None,
				eleKey=None, structKey=None, charge=None):
		""" Initializer for class used to create the standardinput object for grid-convergence calculations
		
		Args:
			geom: (plato_pylib UnitCell object) geometry to use
			kPts: (len-3 int iter) Monkhorst-Pack k-points to use. e.g [20,14,6] means 20 along x, 14 along y, 6 along z
			maxScf: (int) Number of scf cycles to use. One is the tutorial recommended value
			cp2kMethodStr: (str) String used to request a basic pyCP2K object. See cp2k.method_register.getRegisteredCP2KObjCreatorStrings() for options (or to add your own)
			basisKeyDicts: (iter of dicts) This contains all neccesary information on the basis sets we need. Each entry is a dictionary of the form {eleKey:basisKey}. Each eleKey is an elements chemical symbol (e.g. Mg) . Each basisKey corresponds to a key used to obtain a basis set from cp2k.basis_register
			basisAlias: (iter of strs) Each entry is an alias to attach to a basis set defined in basisKeyDicts. This is used to differentiate between them in the output and used to keep folders separate
			workFolder: (str) Path to a folder to run the calculations in. Various sub-folders will be created from here to keep the calculations separate
			eleKey: (str) Simply used to label the element(or chemical compound) we're looking at. Must be set, but its value likely wont matter
			structKey: (str) Just used to label the structure. If your only looking at 1 structure (likely) then it doesnt really matter what you set this to (but it MUST be set)
			workFlowToOutputObjMap: (function, f(workflow)). Takes the gridConvergence workflow and extract data you want in the format you want. If this is left blank then standardOutput.data will simply be workflow.output. Setting to the mapWorkflowOutputToArray function means the output data will be an nx2 array of convergence values vs delta Energy values (relative to energy of the most converged case)
			charge: (int) The charge on the system (default=0)

		Note:
			The createStandardInput function is whats actually used to get a StandardInput object

		"""

		self.geom = geom
		self.kPts = kPts
		self.maxScf = maxScf
		self.methodLabel = cp2kMethodStr
		self.basisKeyDicts = basisKeyDicts #iter is best
		self.basisAliases = basisAliases #iter is best
		self.baseWorkfolder = workFolder
		self.eleKey = eleKey
		self.structKey = structKey
		self.charge = charge
	
		self.workFlowToOutputObjMap = workFlowToOutputObjMap
	
		self._gridValsToUse = None
	
	def createStandardInput(self, gridValsToUse):
		self._gridValsToUse = gridValsToUse
		outObjs = list()
		for basisKey, basisAlias in it.zip_longest(self.basisKeyDicts, self.basisAliases):
			outObjs.append( self._createStandardInputForOneBasis(basisKey,basisAlias) )
		return calcRunners.StandardInputObjComposite(outObjs)
		
	def _createStandardInputForOneBasis(self, basisKeyDict, basisAlias):
		#Figure out the type of grid-convergence we're using
		gridVaryVals = self._getGridVaryVals()
		gridConvType = self._getGridConvType()

		#Create ALL the CP2K objects we need 
		calcFolder = os.path.join(self.baseWorkfolder, basisAlias, gridConvType)
		outLabel = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=basisAlias)
		allCP2KObjs = self._getAllRequiredCP2KObjs( basisKeyDict, calcFolder)

		#Create the input object
		outWorkflow = wflowConv.GridConvergenceEnergyWorkflow(allCP2KObjs, gridVaryVals)
		inputObj = calcRunners.StandardInputObj(outWorkflow, outLabel, mapFunction=self.workFlowToOutputObjMap)

		return inputObj

	def _getGridVaryVals(self):
		gridConvTypes = [x.varyType for x in self._gridValsToUse]
		assert all([x==gridConvTypes[0] for x in gridConvTypes]), "Grid-vals should all be for the same convergence type"
		gridConvType = gridConvTypes[0]
		gridVaryVals = [x.varyValCut for x in self._gridValsToUse]
		return gridVaryVals

	def _getGridConvType(self):
		gridConvTypes = [x.varyType for x in self._gridValsToUse]
		assert all([x==gridConvTypes[0] for x in gridConvTypes]), "Grid-vals should all be for the same convergence type"
		gridConvType = gridConvTypes[0]
		return gridConvType


	def _createBaseCP2KObj(self, basisKeyDict):
		basicObj = self._createBasicObject()
		basisObjects = self._getBasisObjects(basisKeyDict)

		pyCP2KHelp.modCp2kObjBasedOnDict(basicObj, {"kpts":self.kPts, "charge":self.charge} )
		pyCP2KHelp.addGeomAndBasisInfoToSimpleCP2KObj(basicObj, self.geom, basisObjects)
		baseCP2KObj = baseCalcObjCP2K.CP2KCalcObj( basicObj )
		baseCP2KObj.maxScf = self.maxScf
		return baseCP2KObj	


	
	def _getBasisObjects(self, basisKeyDict):
		outList = list()
		for eleKey, basisKey in basisKeyDict.items():
			currBasisObj = basRegister.createCP2KBasisObjFromEleAndBasisStr(eleKey,basisKey)
			outList.append( currBasisObj )
		return outList
	
	def _getAllRequiredCP2KObjs(self, basisKeyDict, startFolder):
		allObjs = list()
		for gVals in self._gridValsToUse:
			currObj = self._createBaseCP2KObj(basisKeyDict) 
			currObj.absGridCutoff = gVals.absCut
			currObj.relGridCutoff = gVals.relCut
			currObj.basePath = os.path.join(startFolder, gVals.toStr() )
			allObjs.append(currObj)
		return allObjs
	
	#Not including geom i guess. Need to add geom and k-points separately
	def _createBasicObject(self):
		baseObj = methRegister.createCP2KObjFromMethodStr(self.methodLabel)
		return baseObj


	

def mapWorkflowOutputToArray(inpObj):
	inpObj.workflow.run() #Need to do this to populate the output object
	assert len(inpObj.workflow.namespaceAttrs)==1, "workflow.namespaceAttrs = {}, but there should only be a single output"
	outAttr = inpObj.workflow.namespaceAttrs[0]
	outVals = getattr(inpObj.workflow.output,outAttr) #This should be and [(x,y)] object, i.e. list of tuples with data
	
	#Now we put it in delta eV, relative to energy of the most converged case
	maxIndex, unused = max(enumerate([x[0] for x in outVals]), key=lambda x:x[1])
	refEnergy = outVals[maxIndex][1]
	outVals = [(x,y-refEnergy) for x,y in outVals]
	
	return outVals


