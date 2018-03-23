from enum import IntEnum


class DepictionSource(IntEnum):
    Pubchem = 1
    Template = 2
    RdKit = 3


class ConformerType(IntEnum):
    Ideal = 1
    Model = 2
    Depiction = 3
    Computed = 4
