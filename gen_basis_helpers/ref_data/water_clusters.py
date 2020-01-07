
import itertools as it
from . import base_objs as baseObjs
from . import helpers_ref_data as helpers
from ..shared import unit_convs as uConv
from ..shared import label_objs as labelObjs
from ..shared import misc_utils as misc

ALL_WATER_CLUSTER_KEYS = set()
WATER_CLUSTERS_FUNCTION_GETTER_DICT = dict()

#Key should be a tuple (geomMethod, energiesMethod, numbWater)
def _registerKeyToWaterClusterDict(key):
	def decorate(funct):
		WATER_CLUSTERS_FUNCTION_GETTER_DICT[key] = funct
		ALL_WATER_CLUSTER_KEYS.add(key)
		return funct
	return decorate


class WaterClusterData(baseObjs.MolecularClusterData):

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=False)
	def __init__(self, label, geom, nMolecules, bindingEnergy):
		self._label  = label
		self._geom = geom
		self._nMolecules = nMolecules
		self._bindingEnergy = bindingEnergy

	@property
	def label(self):
		return [self._label]

	@property
	def geom(self):
		return [self._geom]

	@property
	def nMolecules(self):
		return [self._nMolecules]

	@property
	def bindingEnergy(self):
		return [self._bindingEnergy]


def createWaterClusterRefData(geom="ccsd", energies="pbe-plane-wave", types=["all"], numbWater=6):
	""" Factory function to get geometries and binding energies for water clusters
	
	Args:
		geom (str): The method used to optimise cluster geometries [opts = ccsd]
		energies (str): The method used to calculate binding energies [opts = pbe-plane-wave, None]
		types (list): Names of structrues you want included. If "all" is present in the list then all possible will be included
		numbWater(int): Number of water molecules in the clusters requested

	Returns
		 CompositeMolecularClusterData object containing geometries/labels/binding energies for all requested
 
	Raises:
		 KeyError: If any of the options (geom/energies/types/numbWater) are inconsistent (missing)
	"""

	if "all" in types:
		geomStrs = _getAllWaterClusterTypeStrs(numbWater)
	else:
		geomStrs = types

	comboKey = (geom, energies, numbWater)

	#Get all the geometries, then filter out all those not in geomStrs
	baseComposite = WATER_CLUSTERS_FUNCTION_GETTER_DICT[comboKey]()
	relevantObjs = list()
	for x in geomStrs:
		currObjs = baseComposite.getObjectsWithComponents([x]) #Should be a length 1 list really
		print("x = {}, currObjs = {}".format(x,currObjs))
		relevantObjs.extend(currObjs)

	print("The value of relevantObjs = {}".format(relevantObjs))
	outComposite = baseObjs.CompositeMolecularClusterData( relevantObjs )

	return outComposite


def _getAllWaterClusterTypeStrs(numbWater:int):
	if numbWater == 1:
		return ["monomer"]
	elif numbWater == 6:
		return ["prism","book", "ring", "cage"]
	else:
		raise ValueError("No types of clusters known for numbWater = {}".format(numbWater))

	return None


@_registerKeyToWaterClusterDict( ("ccsd", "pbe-plane-wave", 1) )
def _createWaterClusterDataMonomer_ccsdGeoms_planeWaveEnergies():
	structName = "monomer"
	energy = 0.0
	geom = _getMonomerGeom_ccsd()
	nMolecules = 1
	label = labelObjs.StandardLabel(eleKey="h20", methodKey="plane-wave", structKey="monomer")

	return baseObjs.CompositeMolecularClusterData( [WaterClusterData(label, geom, nMolecules, energy)] )

@_registerKeyToWaterClusterDict( ("ccsd", "pbe-plane-wave", 6) )
def _createWaterClusterDataHexamers_ccsdGeoms_planeWaveEnergies():
	structNames = ["prism","book", "ring", "cage"]
	energies = [0.0,0.0,0.0,0.0] #TODO: Need to actually calculate these
	geoms = [_getHexamerPrismGeom_ccsd(), _getHexamerBookGeom_ccsd(), _getHexamerRingGeom_ccsd(), getHexamerCageGeom_ccsd()]
	nMolecules = [6 for x in structNames]
	labels = [labelObjs.StandardLabel(eleKey="h2o6", methodKey="plane-wave", structKey=struct) for struct in structNames]

	outObjs = list()
	for label, geom, energy in it.zip_longest(labels,geoms,energies):
		outObjs.append( WaterClusterData(label,geom,6, energy) )

	return baseObjs.CompositeMolecularClusterData(outObjs)


def _getMonomerGeom_ccsd():
	""" Taken from 10.1063/1.4922262 ESI (CCSD(T) aug-cc-pV5Z basis) """

	boxSize = 30
	cartCoords = [ [0.00000, 0.00000, 0.95843],
	               [0.00000, 0.00000, 0.00000],
	               [0.00000, 0.92823, 1.19715] ]

	atomSymbols = ["O","H","H"]

	return _getWaterClusterUCellGeomFromCartCoordsAndSymbolList( cartCoords, atomSymbols, boxSize=boxSize)

def _getHexamerBookGeom_ccsd():
	""" Taken from 10.1063/1.4922262 ESI (CCSD(T) aug-cc-pVDZ basis) """

	boxSize = 30
	#These are in angstrom
	cartCoords = [ [1.56069, -0.12812, 0.88613],
	               [1.52699, -0.97048, 0.3703],
	               [2.1873, -0.28764, 1.60295],
	               [1.20059, -2.43241, -0.49624],
	               [0.21531, -2.49121, -0.52049],
	               [1.47832, -2.52835, -1.41595],
	               [-1.51245, -2.26728, -0.45096],
	               [-2.08299, -2.88888, 0.0178],
	               [-1.58672, -1.41889, 0.05059],
	               [-1.4059, 2.47953, -0.43913],
	               [-0.43051, 2.5434, -0.51185],
	               [-1.67599, 3.29656, -0.00125],
	               [1.38351, 2.315, -0.46487],
	               [1.55876, 1.45168, -0.03765],
	               [1.88494, 2.29473, -1.2896],
	               [-1.3686, 0.03998, 0.97301],
	               [-0.41494, 0.05962, 1.16481],
	               [-1.52641, 0.88578, 0.50611] ]

	_convertCartCoordsFromAngToBohr(cartCoords)
	atomSymbols = ["O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H"]

	return _getWaterClusterUCellGeomFromCartCoordsAndSymbolList( cartCoords, atomSymbols, boxSize=boxSize)


def _getHexamerPrismGeom_ccsd():
	""" Taken from 10.1063/1.4922262 ESI (CCSD(T) aug-cc-pVDZ basis) """

	cartCoords = [ [-1.22193, -1.90531, -0.21831],
	               [-1.94971, -2.52152, -0.3679],
	               [-1.6294, -1.00895, -0.22518],
	               [1.04638, -1.13295, 1.45658],
	               [1.42898, -1.13757, 0.55267],
	               [0.22407, -1.6311, 1.30815],
	               [0.26034, 1.42866, 1.49801],
	               [0.59137, 0.49653, 1.60457],
	               [0.47875, 1.88448, 2.32003],
	               [1.36252, -1.04491, -1.33461],
	               [0.48734, -1.47004, -1.31362],
	               [1.15792, -0.10614, -1.49622],
	               [0.58342, 1.77029, -1.24951],
	               [1.08164, 2.50916, -1.62138],
	               [0.65773, 1.86277, -0.27756],
	               [-1.97606, 0.74641, -0.17154],
	               [-1.51212, 1.15948, -0.91783],
	               [-1.45544, 1.06499, 0.5888] ]

	_convertCartCoordsFromAngToBohr(cartCoords)
	atomSymbols = ["O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H"]

	return _getWaterClusterUCellGeomFromCartCoordsAndSymbolList( cartCoords, atomSymbols)

def _getHexamerRingGeom_ccsd():
	""" Taken from 10.1063/1.4922262 ESI (CCSD(T) aug-cc-pVDZ basis) """

	cartCoords = [ [2.72357, -7E-05, -0.1395],
	               [2.22243, -0.83964, -0.02832],
	               [3.23444, -0.12519, -0.94897],
	               [-1.36172, 2.35871, -0.1395],
	               [-0.38406, 2.3445, -0.02832],
	               [-1.5088, 2.8637, -0.94897],
	               [-1.36184, -2.35864, -0.1395],
	               [-1.83837, -1.50486, -0.02832],
	               [-1.72564, -2.73851, -0.94897],
	               [-2.72357, 7E-05, 0.1395],
	               [-2.22243, 0.83964, 0.02832],
	               [-3.23444, 0.12519, 0.94897],
	               [1.36172, -2.35871, 0.1395],
	               [0.38406, -2.3445, 0.02832],
	               [1.5088, -2.8637, 0.94897],
	               [1.36184, 2.35864, 0.1395],
	               [1.83837, 1.50486, 0.02832],
	               [1.72564, 2.73851, 0.94897] ]
	
	_convertCartCoordsFromAngToBohr(cartCoords)
	atomSymbols = ["O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H"]

	return _getWaterClusterUCellGeomFromCartCoordsAndSymbolList( cartCoords, atomSymbols)


def getHexamerCageGeom_ccsd():
	""" Taken from 10.1063/1.4922262 ESI (CCSD(T) aug-cc-pVDZ basis) """

	cartCoords = [ [-1.661, -0.88941, -0.61874],
	               [-1.16356, -1.70986, -0.40113],
	               [-2.54803, -1.17871, -0.86557],
	               [-0.7476, 0.79844, 1.59127],
	               [-1.2538, 0.24714, 0.96714],
	               [0.11957, 0.35685, 1.60517],
	               [0.60262, 0.66269, -1.60343],
	               [-0.25935, 0.21177, -1.5587],
	               [0.40689, 1.56097, -1.26686],
	               [1.66625, -0.56398, 0.55808],
	               [2.61024, -0.41235, 0.68959],
	               [1.42292, -0.07583, -0.27131],
	               [0.11793, -2.82571, 0.17266],
	               [0.51806, -3.46631, -0.42902],
	               [0.8239, -2.16772, 0.34894],
	               [0.04554, 2.9014, -0.02762],
	               [-0.3343, 2.2942, 0.64583],
	               [-0.53282, 3.67327, -0.04122] ]

	_convertCartCoordsFromAngToBohr(cartCoords)
	atomSymbols =  ["O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H", "O", "H", "H"] 

	return _getWaterClusterUCellGeomFromCartCoordsAndSymbolList( cartCoords, atomSymbols)


def _getWaterClusterUCellGeomFromCartCoordsAndSymbolList(cartCoords, atomSymbols, boxSize=30):
	""" Cart-coords should be a length n iter of 3-element lists, atomSymbols should be a length n iter of atomic symbols """
	fullCart = [coords + [symbol] for coords,symbol in it.zip_longest(cartCoords, atomSymbols)]
	outUCell = helpers.getEmptyCubicBoxUCell(boxSize)
	outUCell.cartCoords = fullCart
	return outUCell

def _convertCartCoordsFromAngToBohr(cartCoords):
	""" Cart-Coords should be an iter of 3-element lists (x,y,z) """
	for idx,coords in enumerate(cartCoords):
		cartCoords[idx] = [x*uConv.ANG_TO_BOHR for x in coords]

