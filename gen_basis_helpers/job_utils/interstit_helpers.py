
import plato_pylib.utils.supercell as supCell
import plato_fit_integrals.core.workflow_coordinator as wFlowCoordinator
import plato_fit_integrals.initialise.create_interstit_workflows as createWFlows
import plato_pylib.utils.job_running_functs as jobRun


class InterstitialType():
	""" Defines a type of interstitial independent of element and method used. """
	def __init__(self, relaxed, structType, interstitType, cellDims, runCalcs):
		""" Description of function
		
		Args:
			relaxed: str, unrelaxed/relaxed_constant_p/relaxed_constant_v
			structType: str, the base crystal type (hcp/bcc/fcc)
			interstitType: str describing the deformation (e.g. octahedral)
			cellDims: iter(3 element), number of unit cells in each direction
			runCalcs: Bool, Whether to run calculations for this 
		"""

		self.relaxed = relaxed
		self.structType = structType
		self.interstitType = interstitType
		self.cellDims = cellDims
		self.runCalcs = runCalcs
		
	#Could maybe pull this stuff to its onw inter_helpers to make testing simpler
	def getInterStructFromRefDataStruct(self,refDataStruct):
		""" Gets the relevant structures using a reference object (
		
		Args:
			refDataStruct: RefElementalDataBase object for one element
				
		Returns
			outStructs: UCell structures 
		
		Raises:
			Errors
		"""

		cellStr = "_".join([str(x) for x in self.cellDims])
		outObj = refDataStruct.getSelfInterstitialPlaneWaveStruct(self.structType, self.interstitType, self.relaxed, cellStr)
		return outObj





class InterstitialCalculationsComposite():
	def __init__(self, inpObjs, label):
		self._objs = inpObjs
		self.label = label
		

	def runCalcs(self):
		jobRun.executeRunCommsParralel(self.runComms,self.nCores,quiet=False)

	@property
	def nCores(self):
		return max([x.nCores for x in self._objs])

	@property
	def runComms(self):
		allComms = list()
		for x in self._objs:
			allComms.extend(x.runComms)
		return allComms

	def printResults(self):
		print("\n {} \n".format(self.label))
		for x in self._objs:
			x.printResults()


#NOTE: Recommended to control the inclPreRun command at a lower level. Theres an option in the interstitial object creation to 
# not generate run-commands, this can be used to selectivey turn off/on jobs.
class InterstitialCalculationSet():
	def __init__(self, workFlowCoord, label, runJobs=True):
		self.workFlowCoord = workFlowCoord
		self.label = label
		self.runJobs = runJobs
	
	@property
	def nCores(self):
		return self.workFlowCoord.nCores

	@property
	def runComms(self):
		return self.workFlowCoord.preRunShellComms
	
	def runCalcs(self):
		if self.runJobs:
			self.workFlowCoord.run(inclPreRun=True)

	def printResults(self):
		self.workFlowCoord.run(inclPreRun=False)
		print(self.label)
		outVals = self.workFlowCoord.propertyValues
		attrs = vars(outVals)
		[print("{} = {}".format(key, getattr(outVals,key))) for key in attrs]
		print("\n")



class InterstitialCalcSetFactory():
	def __init__(self,workFolder, platoComm, optDict, allInterstitInfoObjs, refData, nCores, label, geomType="plane-wave", eType="electronicCohesiveE"):
		""" Initialise factory instance.
		
		Args:
			workFolder: str, path to base folder to carry out calculations in
			platoComm: str, represents plato program to run (tb1/dft2/dft)
			optDict: dict, plato_pylib style input options dict for plato (only needs to contain options you want changed from default)
			allInterstitInfoObjs: iter, InterstitialType objects. These contain all info on the interstitials you want calculated (i.e. independent of the method used)
			refData: RefElementalDataBase object for one element
			nCores: Number of cores to use
			label: Str used to differentiate this object from others. Example might be the method used (e.g. tb1_2c)
			geomType(Optional): str, what type of geometry to use for the reference structure. currently only plane-wave (will add expt later)
			eType(Optional): str, the energy type to use. See energies object in plato_pylib. One possible option is electronicTotalE
		"""

		self.workFolder = workFolder
		self.platoComm = platoComm
		self.optDict = optDict
		self.interObjs = list(allInterstitInfoObjs)
		self.refData = refData	
		self.geomType = geomType
		self.nCores = nCores
		self.label = label
		self.eType = eType

	def _getRefGeom(self, interObj):
		if self.geomType == "plane-wave":
			return self.refData.getPlaneWaveGeom(interObj.structType)
		elif self.geomType == "expt":
			return self.refData.getExptGeom(interObj.structType)
		else:
			raise ValueError("{} is an invalid options for geomType".format(self.geomType))


	def _createWorkFlowCoordinator(self):
		allObjs = list()
		for x in self.interObjs:
			cellDims = x.cellDims
			interStruct = x.getInterStructFromRefDataStruct(self.refData)
			refStruct = self._getRefGeom(x) #Should be the same for all really
			refStruct = supCell.superCellFromUCell(refStruct, cellDims)
			startFolder = self.workFolder
			modOptDict = self.optDict
			platoComm = self.platoComm
			relaxed = x.relaxed
			interType = x.interstitType
			runCalcs = x.runCalcs
			currObj = createWFlows.CreateInterstitialWorkFlow(refStruct, interStruct, startFolder, modOptDict, platoComm, relaxed=relaxed, cellDims=cellDims, interType=interType, genPreShellComms=runCalcs, eType=self.eType)()
			allObjs.append(currObj)

		workFlowCoord = wFlowCoordinator.WorkFlowCoordinator(allObjs,nCores=self.nCores,quietPreShellComms=False)
		return workFlowCoord

	def __call__(self):
		wFlowCoord = self._createWorkFlowCoordinator()
		outObj = InterstitialCalculationSet(wFlowCoord, self.label)
		return outObj




