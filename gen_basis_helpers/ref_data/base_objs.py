
from ..shared import misc_utils as misc


class MolecularClusterData():

	@property
	def geom(self):
		""" UnitCell object containing geometry of cluster (list returned to enable composite pattern)
		"""
		raise NotImplementedError()

	@property
	def nMolecules(self):
		""" Number of molecules in the cluster (integer list returned to enable composite pattern)
		"""
		raise NotImplementedError()

	@property
	def bindingEnergy(self):
		""" Total binding energy of the cluster ( E(cluster) - E(monomer)*nMoleculesInCluster ). List returned to enable composite pattern
		"""
		raise NotImplementedError()

	@property
	def label(self):
		""" StandardLabel object used to identify the cluster this belongs to. List returned to enable composite pattern
		"""
		raise NotImplementedError()



class CompositeMolecularClusterData(MolecularClusterData):
	"""Composite of the MolecularClusterData object (i.e. provides uniform access to a group of MolecularClusterDataObjects)

	"""

	nMolecules = misc.StandardComponentDescriptor("nMolecules")
	bindingEnergy = misc.StandardComponentDescriptor("bindingEnergy")
	label = misc.StandardComponentDescriptor("label")
	geom = misc.StandardComponentDescriptor("geom")

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=True)
	def __init__(self, branchObjs):
		self.objs = list(branchObjs) #This attribute is needed for composite search 

