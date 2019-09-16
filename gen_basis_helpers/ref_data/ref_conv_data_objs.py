


class RefConvergenceDatabase():
	""" Class for holding converged grid values/k-points for various structures """

	@property
	def integGridVals(self):
		raise NotImplementedError

	@property
	def kptGridVals(self):
		raise NotImplementedError

class IntegGridConvergencePureElements():
	""" Class for holding integration grid convergence to use for pure elements in various structs """

	def getPrimCellDftGrid(structKey):
		""" Get stored value of grid for primitive cell for dft plato program
		
		Args:
			structKey(str): hcp/bcc/fcc depending on which structure you want
				
		Returns
			dftGrid(float): the FFT spacing value to use in plato dft program for the relevant grid

		Raises:
			KeyError: If the values arent Present

		"""
		raise NotImplementedError

	def getPrimCellDft2AngularGrid(structKey):
		""" Get stored angular grid parameters for primitive cell for dft2 plato program
		
		Args:
			structKey:(str) hcp/bcc/fcc depending on which structure you want
				
		Returns
			gridParams:(len 3 list) contains number of atom-centred points for [radial,theta,phi] 

		Raises:
			KeyError: If the values arent Present

		"""
		raise NotImplementedError


	def getInterstitialDft2AngularGrid(self,structKey, dims):
		""" Get stored integral mesh for an interstitial. Prim-cell values may seg-fault here due to massive atom packing
		
		Args:
			structKey:(str) hcp/bcc/fcc depending on which structure you want
			dims: (3 tuple), dimensions of the supercell in x,y,z	
		Returns
			gridParams:(len 3 list) contains number of atom-centred points for [radial,theta,phi] 

		Raises:

		"""
		raise NotImplementedError

class KPointConvergence():
	""" Class for holding k-points to use for various structs in pure elements """
	

	def getKptsPrimCell(self,structKey):
		""" Get stored k-point parameters for primitive cell 
		
		Args:
			structKey:(str) hcp/bcc/fcc depending on which structure you want
				
		Returns
			kPts:(len 3 list) with number of k-points in each direction for Monkhorst-Pack mesh [kx,ky,kz]

		Raises:
			KeyError: If the values arent Present

		"""
		raise NotImplementedError


	def getKPtsSuperCell(self,structKey,dims):
		""" Get stored k-point parameters for a supercell
		
		Args:
			structKey:(str) hcp/bcc/fcc depending on which structure you want
			dims: (3 tuple), dimensions of the supercell in x,y,z	
		Returns
			kPts:(len 3 list) with number of k-points in each direction for Monkhorst-Pack mesh [kx,ky,kz]

		Raises:
			KeyError: If the values arent Present
		"""
		raise NotImplementedError



