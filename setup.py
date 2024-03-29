from distutils.core import setup

setup(name='gen_basis_helpers',
	  version='1.0',
	  author='Richard Fogarty',
	  author_email = 'richard.m.fogarty@gmail.com',
	  packages = ['gen_basis_helpers', 'gen_basis_helpers.analyse_md', 'gen_basis_helpers.analyse_md.surf_norm_z', 'gen_basis_helpers.analyse_md.specific_papers',
                  'gen_basis_helpers.adsorption' ,'gen_basis_helpers.convergers','gen_basis_helpers.castep','gen_basis_helpers.castep.private',
	              'gen_basis_helpers.cp2k','gen_basis_helpers.cp2k.private', 'gen_basis_helpers.cp2k.ref_sets', 'gen_basis_helpers.cp2k.job_utils',
	              'gen_basis_helpers.create_model', 'gen_basis_helpers.db_help', 'gen_basis_helpers.elemental_eos',
	              'gen_basis_helpers.fit_cp2k_basis', 'gen_basis_helpers.gau_prod_theorem',
	              'gen_basis_helpers.job_utils', 'gen_basis_helpers.job_helpers', 'gen_basis_helpers.job_helpers.castep',
	              'gen_basis_helpers.job_helpers.cp2k','gen_basis_helpers.job_helpers.plato', 
	              'gen_basis_helpers.lammps_interface', 'gen_basis_helpers.misc',
	              'gen_basis_helpers.plato','gen_basis_helpers.ref_data',
	              'gen_basis_helpers.special_builders',
	              'gen_basis_helpers.shared', 'gen_basis_helpers.workflows']
	 )

