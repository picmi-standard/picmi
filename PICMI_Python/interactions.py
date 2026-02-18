"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in this file are related to interactions (e.g. field ionization, collisions, QED)
"""
from __future__ import annotations
from typing import Any

from pydantic import Field

from .base import _ClassWithInit

from .particles import PICMI_Species, PICMI_MultiSpecies

PICMI_SpeciesType = PICMI_Species | PICMI_MultiSpecies

class PICMI_FieldIonization(_ClassWithInit):
    """
    Field ionization on an ion species.
    """
    model: str = Field(description="Ionization model, e.g. 'ADK'")
    ionized_species: "PICMI_Species | PICMI_MultiSpecies" = Field(description="Species that is ionized")
    product_species: "PICMI_Species | PICMI_MultiSpecies" = Field(description="Species in which ionized electrons are stored")

