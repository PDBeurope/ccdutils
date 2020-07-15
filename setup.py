from setuptools import setup, find_namespace_packages
import pdbeccdutils

      install_requires=['Pillow', 'scipy', 'numpy', 'pdbecif'],
setup(
    name="pdbeccdutils",
    version=pdbeccdutils.__version__,
    description="Toolkit to deal with wwPDB chemical components definitions for small molecules.",
    project_urls={
        "Source code": "https://github.com/PDBeurope/ccdutils",
        "Documentation": "https://pdbeurope.github.io/ccdutils/",
    },
    author="Protein Data Bank in Europe",
    author_email="pdbehelp@ebi.ac.uk",
    license="Apache License 2.0.",
    keywords="PDB CCD wwPDB small molecule",
    packages=find_namespace_packages(),
    zip_safe=False,
    include_package_data=True,
    tests_require=["pytest"],
    entry_points={
        "console_scripts": [
            "process_components_cif=pdbeccdutils.scripts.process_components_cif_cli:main",
            "setup_pubchem_library=pdbeccdutils.scripts.setup_pubchem_library_cli:main",
        ]
    },
)
