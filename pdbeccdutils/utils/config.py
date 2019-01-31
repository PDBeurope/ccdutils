"""Package configuration with the paths to the templates
"""
import os
import pdbeccdutils

general_templates = os.path.join(os.path.dirname(pdbeccdutils.__file__), 'data', 'general_templates')
coordgen_templates = os.path.join(os.path.dirname(pdbeccdutils.__file__), 'data', 'coordgen_templates')
fragment_library = os.path.join(os.path.dirname(pdbeccdutils.__file__), 'data', 'fragment_library.tsv')
