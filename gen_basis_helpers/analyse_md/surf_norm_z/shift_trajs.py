

import plato_pylib.shared.ucell_class as uCellHelp


def shiftUnitCellToCentreAverageOfZIndices_trajInterface(trajObj, inpIndices, targZ, foldAfter=True):
	""" Shifts all co-ordinates in a trajectory such that the average z-position of inpIndices is targZ. Original purpose is to remove effect of translation along z when analysing MD
	
	Args:
		trajObj: (TrajectoryBase object)
		inpIndices: (iter of ints) The indices of the atoms we want centred
		targZ: (float) Target z-value 
		foldAfter: (Bool) If True the final co-ords are folded back into the unit cell

	Returns
		 Nothing; works in place
 
	"""
	for tStep in trajObj:
		_shiftUnitCellToCentreAverageZOfIndices(tStep.unitCell, inpIndices, targZ, foldAfter=foldAfter)


def _shiftUnitCellToCentreAverageZOfIndices(inpCell, inpIndices, targZ, foldAfter=True):
	""" Shifts cartesian co-ords such that the average z-position of inpIndices is targZ. Original purpose is to remove effect of translation along z when analysing MD
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		inpIndices: (iter of ints) The indices of the atoms we want centred
		targZ: (float) Target z-value 
		foldAfter: (Bool) If True the final co-ords are folded back into the unit cell
 
	Returns
		Nothing; works in place
 
	"""

	#Figure out the shift value
	cartCoords = inpCell.cartCoords
	targZVals = [currCoord[2] for idx,currCoord in enumerate(cartCoords) if idx in inpIndices]
	avInputZVal = sum(targZVals) / len(targZVals)
	shiftZVal = targZ - avInputZVal

	#Apply the translation
	tVector = [0,0,shiftZVal]
	uCellHelp.applyTranslationVectorToCartCoords(inpCell, tVector, foldInAfter=foldAfter)




