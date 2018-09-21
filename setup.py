from setuptools import setup

setup(name='pdbeccdutils',
      version='0.2',
      description='Toolkit to deal with wwPDB chemical components definitions for small molecules.',
      project_urls={
          'Source code': 'https://gitlab.ebi.ac.uk/pdbe/ccdutils',
          'Documentation': 'https://pdbe.gitdocs.ebi.ac.uk/ccdutils/',
      },
      author='Protein Data Bank in Europe',
      author_email='pdbehelp@ebi.ac.uk',
      license='Apache License 2.0.',
      keywords='PDB CCD wwPDB small molecule',
      packages=['pdbeccdutils'],
      scripts=['bin/setup_pubchem_library',
               'bin/process_components_cif',
               'bin/protein_interactions'],
      zip_safe=False,
      include_package_data=True,
      install_requires=['Pillow', 'scipy', 'sphinx', 'sphinx_rtd_theme', 'pdbecif'],
      tests_require=['pytest']
      )
