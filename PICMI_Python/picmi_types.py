"""Types definitions for the picmi standard
"""
import typing
from collections.abc import Sequence

import numpy as np

VectorFloat1 = typing.NewType('VectorFloat1', typing.Union[Sequence[float], np.ndarray[1, np.float64]])
VectorInt1 = typing.NewType('VectorInt1', typing.Union[Sequence[int], np.ndarray[1, np.int64]])
VectorString1 = typing.NewType('VectorString1', Sequence[str])

VectorFloat2 = typing.NewType('VectorFloat2', typing.Union[Sequence[float], np.ndarray[2, np.float64]])
VectorInt2 = typing.NewType('VectorInt2', typing.Union[Sequence[int], np.ndarray[2, np.int64]])
VectorString2 = typing.NewType('VectorString2', Sequence[str])

VectorFloat3 = typing.NewType('VectorFloat3', typing.Union[Sequence[float], np.ndarray[3, np.float64]])
VectorFloatNone3 = typing.NewType('VectorFloatNone3', typing.Union[Sequence[typing.Union[float, None]], np.ndarray[3, np.float64]])
VectorInt3 = typing.NewType('VectorInt3', typing.Union[Sequence[int], np.ndarray[3, np.int64]])
VectorBool3 = typing.NewType('VectorBool3', Sequence[bool])
VectorString3 = typing.NewType('VectorString3', Sequence[str])
Expression = typing.NewType('Expression', str)
VectorExpression3 = typing.NewType('VectorExpression3', Sequence[str])


# These must be defined as strings to avoid circular definitions
GridType = """typing.NewType('GridType', typing.Union[fields.PICMI_Cartesian1DGrid,
                                                      fields.PICMI_CylindricalGrid,
                                                      fields.PICMI_Cartesian2DGrid,
                                                      fields.PICMI_Cartesian3DGrid])"""
# This is needed in the fields.py file since it cannot import itself
FieldsGridType = """typing.NewType('GridType', typing.Union[PICMI_Cartesian1DGrid,
                                                            PICMI_CylindricalGrid,
                                                            PICMI_Cartesian2DGrid,
                                                            PICMI_Cartesian3DGrid])"""

SmootherType = """typing.NewType('SmootherType', PICMI_BinomialSmoother)"""

_species_union = "typing.Union[particles.PICMI_Species, particles.PICMI_MultiSpecies]"

SpeciesType = f"""typing.NewType('SpeciesType', {_species_union})"""

# Note that in the strings, types defined in this file cannot be used, so the _species_union needs to
# be included explicitly.
# The SpeciesArgument allows either a single species or a sequence of species.
SpeciesArgument = f"""typing.NewType('SpeciesArgument', typing.Union[{_species_union},
                                                                     Sequence[{_species_union}]])"""

# Note that this is only used in particles.py and uses classes defined there.
# The "particles" prefix can not be used on the classes since the module
# should not reference itself.
DistributionType = """typing.NewType('DistributionType', typing.Union[PICMI_GaussianBunchDistribution,
                                                                      PICMI_UniformDistribution,
                                                                      PICMI_AnalyticDistribution,
                                                                      PICMI_ParticleListDistribution])"""

LayoutType = """typing.NewType('LayoutType', typing.Union[particles.PICMI_GriddedLayout,
                                                          particles.PICMI_PseudoRandomLayout])"""

SolverType = """typing.NewType('SolverType', typing.Union[fields.PICMI_ElectromagneticSolver,
                                                          fields.PICMI_ElectrostaticSolver,
                                                          fields.PICMI_MagnetostaticSolver])"""

LaserType = """typing.NewType('LaserType', typing.Union[lasers.PICMI_GaussianLaser,
                                                        lasers.PICMI_AnalyticLaser])"""

LaserInjectionType = """typing.NewType('LaserInjectionType', lasers.PICMI_LaserAntenna)"""

AppliedFieldType = """typing.NewType('AppliedFieldType', typing.Union[applied_fields.PICMI_ConstantAppliedField,
                                                                      applied_fields.PICMI_AnalyticAppliedField])"""

DiagnosticType = """typing.NewType('DiagnosticType', typing.Union[diagnostics.PICMI_FieldDiagnostic,
                                                                  diagnostics.PICMI_ElectrostaticFieldDiagnostic,
                                                                  diagnostics.PICMI_ParticleDiagnostic,
                                                                  diagnostics.PICMI_LabFrameFieldDiagnostic,
                                                                  diagnostics.PICMI_LabFrameParticleDiagnostic])"""

