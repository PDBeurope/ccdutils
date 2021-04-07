import os

from setuptools import find_namespace_packages, setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]

    raise RuntimeError("Unable to find version string.")


setup(
    name="pdbeccdutils",
    version=get_version("pdbeccdutils/__init__.py"),
    description="Toolkit to deal with wwPDB chemical components definitions for small molecules.",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    project_urls={
        "Source code": "https://github.com/PDBeurope/ccdutils",
        "Documentation": "https://pdbeurope.github.io/ccdutils/",
    },
    author="Protein Data Bank in Europe",
    author_email="pdbehelp@ebi.ac.uk",
    license="Apache License 2.0.",
    keywords="PDB CCD wwPDB small molecule",
    url="http://pypi.python.org/pypi/pdbeccdutils/",
    packages=find_namespace_packages(),
    zip_safe=False,
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        "Pillow",
        "scipy",
        "numpy",
        "requests",
        "pdbecif>=1.5",
    ],
    extras_require={
        "tests": ["pytest", "pytest-cov", "pre-commit"],
        "docs": [
            "sphinx",
            "sphinx_rtd_theme",
            "recommonmark",
            "sphinx-autodoc-typehints",
            "sphinx-markdown-tables",
        ],
    },
    entry_points={
        "console_scripts": [
            "process_components_cif=pdbeccdutils.scripts.process_components_cif_cli:main",
            "setup_pubchem_library=pdbeccdutils.scripts.setup_pubchem_library_cli:main",
        ]
    },
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Development Status :: 5 - Production/Stable",
    ],
)
