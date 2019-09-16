from distutils.core import setup

setup(name='gen_basis_helpers',
	  version='1.0',
	  author='Richard Fogarty',
	  author_email = 'richard.m.fogarty@gmail.com',
	  packages = ['gen_basis_helpers','gen_basis_helpers.cp2k', 'gen_basis_helpers.create_model', 'gen_basis_helpers.gau_prod_theorem', 'gen_basis_helpers.job_utils', 'gen_basis_helpers.ref_data', 'gen_basis_helpers.shared']
	 )

