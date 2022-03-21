"""Types definitions for the picmi standard
"""
import typing
from collections.abc import Sequence

import numpy as np

VectorFloat3 = typing.NewType('VectorFloat3', typing.Union[Sequence[float], np.ndarray[3, np.float64]])
VectorInt3 = typing.NewType('VectorInt3', typing.Union[Sequence[int], np.ndarray[3, np.int64]])
Expression = typing.NewType('Expression', str)


# These must be defined as strings to avoid circular definitions
GridType = """typing.NewType('GridType', typing.Union[fields.PICMI_Cartesian1DGrid,
                                                   fields.PICMI_CylindricalGrid,
                                                   fields.PICMI_Cartesian2DGrid,
                                                   fields.PICMI_Cartesian3DGrid])"""

SpeciesType = """typing.NewType('SpeciesType', typing.Union[particles.PICMI_Species,
                                                         particles.PICMI_MultiSpecies])"""

SpeciesArgument = """typing.NewType('SpeciesArgument', typing.Union[SpeciesType, Sequence[SpeciesType]])"""


