"""Package configuration with the paths to the templates"""

import os
from pathlib import Path


def get_data_dir():
    return Path(os.path.dirname(__file__)).parent / "data"


general_templates = os.path.join(get_data_dir(), "general_templates")
coordgen_templates = os.path.join(get_data_dir(), "coordgen_templates")
fragment_library = os.path.join(get_data_dir(), "fragment_library.tsv")

DISCARDED_RESIDUES = {"HOH", "UNX"}
