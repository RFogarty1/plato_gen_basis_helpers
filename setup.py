from distutils.core import setup

setup(name='gen_basis_helpers',
	  version='1.0',
	  author='Richard Fogarty',
	  author_email = 'richard.m.fogarty@gmail.com',
	  packages = ['gen_basis_helpers','gen_basis_helpers.convergers','gen_basis_helpers.cp2k',
	              'gen_basis_helpers.cp2k.private', 'gen_basis_helpers.cp2k.ref_sets', 'gen_basis_helpers.cp2k.job_utils',
	              'gen_basis_helpers.create_model', 'gen_basis_helpers.elemental_eos','gen_basis_helpers.gau_prod_theorem',
	              'gen_basis_helpers.job_utils', 'gen_basis_helpers.job_helpers', 'gen_basis_helpers.job_helpers.plato', 'gen_basis_helpers.plato','gen_basis_helpers.ref_data', 
	              'gen_basis_helpers.shared', 'gen_basis_helpers.workflows']
	 )

