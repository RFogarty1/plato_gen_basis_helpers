
import types

import plato_pylib.utils.supercell as supCell

from ..shared import label_objs as labelHelp

from ..shared import calc_runners as calcRunners
from ..shared import creator_resetable_kwargs as baseCreator
from ..shared import base_surface as baseSurf
from ..shared import surfaces as surfGetterHelp
from ..shared import label_objs as labelHelp

from ..workflows import surface_energies as surfFlow

class CodeSpecificStandardInputCreatorTemplate(baseCreator.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("baseGeom")
	registeredKwargs.add("cellDims")
	registeredKwargs.add("surfType")
	registeredKwargs.add("nLayers")
	registeredKwargs.add("lenVac")
	registeredKwargs.add("baseWorkFolder")
	registeredKwargs.add("kPts")
	registeredKwargs.add("eleKey")
	registeredKwargs.add("methodKey")
	registeredKwargs.add("structKey")
	registeredKwargs.add("applyNLayersToBulk")


	def _setDefaultInitAttrs(self):
		self.applyNLayersToBulk = False

	def _createFromSelf(self):
		bulkCalcObj = self._getBulkCalcObj()
		surfCalcObj = self._getSurfaceCalcObj()
		surfAreaFromUnitCellFunct = surfFlow.SurfaceAreaFromUnitCellFunct(self._getSurfaceObjClass())
		workflow = surfFlow.SurfaceEnergyWorkflow(surfCalcObj,bulkCalcObj,surfAreaFromUnitCellFunct)
		label = self.label
		return calcRunners.StandardInputObj(workflow,label)

	# This should be overwritten by code-specific versions
	def _createCalcObjCreator(self):
		""" Return a CalcMethodFactoryBase instance with no kPts or geom present
		"""
		raise NotImplementedError("")

	def _getSurfaceCalcObj(self):
		outCreator = self._createCalcObjCreator() #No geometry included
		outCreator.geom = self._surfaceCell
		outCreator.kPts = self._surfKPts
		outCreator.fileName = self._surfFileName
		return outCreator.create()

	def _getBulkCalcObj(self):
		outCreator = self._createCalcObjCreator()
		outCreator.geom = self._bulkCell
		outCreator.kPts = self.kPts
		outCreator.fileName = "bulk_cell_for_n{}".format(self.nLayers)
		return outCreator.create()

	def _getSurfaceObjClass(self):
		if isinstance(self.surfType, str):
			return self._getSurfaceObjClassFromStr(self.surfType)

		#cant use isinstance to compare classes; have to initialiase actual objects while bypassing __init__ to avoid passing the args
		tempObj = self.surfType.__new__(self.surfType)

		if isinstance(tempObj, baseSurf.BaseSurface):
			return self.surfType
		else:
			raise ValueError("{} is an invalid type".format( type(self.surfType) ))


	def _getSurfaceObjClassFromStr(self, inpStr):
		if inpStr == "hcp0001":
			return surfGetterHelp.Hcp0001Surface
		else:
			raise ValueError("{} is an invalid surface type".format(self.surfType))

	@property
	def _surfFileName(self):
		return "surface_n{}_vac_{:.2f}".format(self.nLayers,self.lenVac).replace(".","pt")


	@property
	def _bulkCell(self):
		baseBulkCell = self._bulkCellNoSurfaceLayers
		if self.applyNLayersToBulk:
			surfClass = self._getSurfaceObjClass() #TODO: Need to switch this to a factory soon, to unify interface between diff surfaces
			lenVac = 0.0
			outUCell = surfClass(baseBulkCell,self.nLayers,lenVac).unitCell
		else:
			outUCell = baseBulkCell
		return outUCell

	@property
	def _surfaceCell(self):
		bulkCell = self._bulkCellNoSurfaceLayers
		surfClass = self._getSurfaceObjClass() #TODO: Need to switch this to a factory soon, to unify interface between diff surfaces
		surfUCell = surfClass(bulkCell,self.nLayers,self.lenVac).unitCell
		return surfUCell

	@property
	def _bulkCellNoSurfaceLayers(self):
		return supCell.superCellFromUCell( self.baseGeom, self.cellDims )


	@property
	def _surfKPts(self):
		return [self.kPts[0], self.kPts[1], 1]

	@property
	def label(self):
		return labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)

	#Property only so i can apply a docstring
	@property
	def applyNLayersToBulk(self):
		""" If True then the bulk calculation will simply be on the same system as the surface, except lenVac=0 for bulk. If False we use the minimal cell for the bulk calculation. False is faster, True might sometimes be better if we want to minimize the effects of using different numbers of k-points
		"""
		return self._applyNLayersToBulk

	@applyNLayersToBulk.setter
	def applyNLayersToBulk(self,val):
		self._applyNLayersToBulk = val if val is not None else False



class MapSurfaceEnergiesToStandardFormat():
	"""Callable class used to transform surface energy workflow output into more useful format

	"""

	def __init__(self, ePerAtomFmtStr="{:.3g}", surfEnergyFmtStr="{:.4f}"):
		self.ePerAtomFmtStr = ePerAtomFmtStr
		self.surfEnergyFmtStr = surfEnergyFmtStr

	def _getTableData(self,stdInputObj):
		assert len(stdInputObj.label)==1
		assert len(stdInputObj.workflow.output)==1
		methKey = stdInputObj.label[0].methodKey
		surfEnergy = self.surfEnergyFmtStr.format(stdInputObj.workflow.output[0].surfaceEnergy)
		return [methKey, surfEnergy]

	def _getTableDataWithEPerAtom(self,stdInputObj):
		outTable = self._getTableData(stdInputObj)
		ePerAtomBulk = stdInputObj.workflow.output[0].bulkEPerAtom
		ePerAtomSurf = stdInputObj.workflow.output[0].surfEPerAtom
		outTable += [self.ePerAtomFmtStr.format(x) for x in [ePerAtomBulk,ePerAtomSurf]]  
		return outTable	


	def _getTableHeadings(self):
		return ["Basis Set", "Surface Energy $eV a_{0}^{-2}$"]

	def _getTableWithEPerAtomHeadings(self):
		return ["Basis Set", "Surface Energy $eV a_{0}^{-2}$", "E per atom (bulk, eV)", "E per atom (surface, eV)",]

	def __call__(self, stdInputObj):
		stdInputObj.workflow.run()
		output = types.SimpleNamespace(tableData=None)
		output.tableData = self._getTableData(stdInputObj)
		output.tableHeaders = self._getTableHeadings()
		output.tableWithEPerAtomVals = self._getTableDataWithEPerAtom(stdInputObj)
		output.tableHeadersWithEPerAtom = self._getTableWithEPerAtomHeadings()
		return output
