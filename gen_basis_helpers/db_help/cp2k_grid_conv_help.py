
import itertools as it
import functools

import plato_pylib.shared.ucell_class as uCellHelp

def getGridConvVsTotalEnergy(varyAttr, queryDict, collection, deltaE=True, ePerAtom=True, eType="electronicTotalE"):
	""" Gets a list of grid convergence values vs total energies
	
	Args:
		varyAttr: (str) Used to get x-values from a record, usually "absGridConv" or "relGridConv"
		queryDict: (dict) Key,vals used to query the relevant database, basically passed to collection.find
		collection: Pymongo collection object
		deltaE: (bool) If True then values are returned relative to that at the HIGEST x-value
		ePerAtom: (bool) If True then energies per atom are returned
		eType: (str) The energy type to get from records["energies"] (which should be a dict representing the plato_pylib energies class
			
	Returns
		outVals: iter of [gridCutoff,energy] (Will be sorted by gridcutoff in ascending order) 
	
	Raises:
		AssertionError: Should be raised if we find no records
	"""
	recordsToOutVals = functools.partial(getEnergiesFromRecords, ePerAtom=ePerAtom, eType="electronicTotalE")
	postProcessFuncts = [_sortOutputByXVals]
	if deltaE:
		postProcessFuncts.append(_getYRelativeToHighestXVal)
	convFunct = GetCP2KGridConvVsProp(varyAttr, recordsToOutVals, postProcessFuncts)
	return convFunct(queryDict, collection)

class GetCP2KGridConvVsProp():

	def __init__(self, varyAttr, recordsToOutVals, postProcessFuncts):
		""" Initializer
		
		Args:
			varyAttr: (str) absGridCutoff or relGridCutoff generally
			recordsToOutVals: f(records)->vals Where records and vals are equal length iters
			postProcessFuncts: iter of 	f([[x,y] for x,y in records])->[[x,y] for x,y in records] functions. Expected to do things like ensure output is sorted or relative to a specific x-val
		"""
		self.varyAttr = varyAttr
		self.recordsToOutVals = recordsToOutVals
		self.postProcessFuncts = list(postProcessFuncts)

	def getVals(self, queryDict, collection):
		#1) Get basic x vs y
		records = [x for x in collection.find(queryDict)]
		assert len(records)>0, "No records found for queryDict = {}".format(queryDict)
		xVals = [x[self.varyAttr] for x in records]
		yVals = self.recordsToOutVals(records)
		outVals = [ [x,y] for x,y in it.zip_longest(xVals, yVals) ]

		#2) Apply any post processing required
		for pFunct in self.postProcessFuncts:
			outVals = pFunct(outVals)

		return outVals

	def __call__(self, queryDict, collection):
		return self.getVals(queryDict, collection)


def getEnergiesFromRecords(records, ePerAtom=False, eType="electronicTotalE"):
	outEnergies = list()
	for record in records:
		currEnergy = record["energies"][eType]
		if ePerAtom:
			nAtoms = len( uCellHelp.UnitCell.fromDict(record["out_geom"]).fractCoords )
			currEnergy = currEnergy / nAtoms
		outEnergies.append(currEnergy)
	return outEnergies

def _sortOutputByXVals(xVsYVals):
	return sorted( xVsYVals, key=lambda x:x[0] )

def _getYRelativeToHighestXVal(xVsYVal):
	labelledXVals = [ [idx,val[0]] for idx,val in enumerate(xVsYVal) ]
	refIdx = max(labelledXVals,key=lambda x:x[1])[0] 
	refVal = xVsYVal[refIdx][1]
	outVals = list()
	for (x,y) in xVsYVal:
		outVals.append( [x, y-refVal] )
	return outVals






