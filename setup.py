from setuptools import setup

setup(name='pdbeccdutils',
      version='0.1.0',
      description='Tools to deal with wwPDB chemical components definitions for small molecules',
      url='https://gitlab.com/pdbe/ccd_utils',
      author='Protein Data Bank in Europe',
      author_email='pdbehelp@ebi.ac.uk',
      license='Apache License 2.0.',
      packages=['pdbeccdutils'],
      scripts=['bin/process_components_cif',
               'bin/ccd_utils_cli'],
      zip_safe=False,
      include_package_data=True,
      install_requires=['yattag', 'lxml', 'pandas']
      )
