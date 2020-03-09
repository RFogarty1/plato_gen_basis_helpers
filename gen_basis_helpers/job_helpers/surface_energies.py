

import plato_pylib.utils.supercell as supCell

from ..shared import label_objs as labelHelp

from ..shared import calc_runners as calcRunners
from ..shared import creator_resetable_kwargs as baseCreator
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
		bulkCalcObj = self._getBulKCalcObj()
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
		if self.surfType == "hcp0001":
			return surfGetterHelp.Hcp0001Surface
		else:
			raise ValueError("{} is an invalid surface type".format(self.surfType))

	@property
	def _surfFileName(self):
		return "surface_n{}_vac_{:.2f}".format(self.nLayers,self.lenVac).replace(".","pt")


	@property
	def _bulkCell(self):
		baseBulkCell = supCell.superCellFromUCell( self.baseGeom, self.cellDims )
		if self.applyNLayersToBulk:
			surfClass = self._getSurfaceObjClass() #TODO: Need to switch this to a factory soon, to unify interface between diff surfaces
			lenVac = 0.0
			outUCell = surfClass(baseBulkCell,self.nLayers,lenVac).unitCell
		else:
			outUCell = baseBulkCell
		return outUCell

	@property
	def _surfaceCell(self):
		bulkCell = self._bulkCell
		surfClass = self._getSurfaceObjClass() #TODO: Need to switch this to a factory soon, to unify interface between diff surfaces
		surfUCell = surfClass(bulkCell,self.nLayers,self.lenVac).unitCell
		return surfUCell

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

