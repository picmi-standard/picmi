"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in this file are related to interactions (e.g. field ionization, collisions, QED)
"""

from pydantic import Field

from .base import PICMI_Interaction, _PICMIModel


class PICMI_FieldIonization(PICMI_Interaction):
    """
    Field ionization on an ion species
    """
    model: str = Field(
        description='Ionization model, e.g. "ADK"'
    )
    # The species class is defined in particles.py, which rebuilds this class once it is defined.
    ionized_species: "PICMI_Species" = Field(
        description="Species that is ionized"
    )
    product_species: "PICMI_Species" = Field(
        description="Species in which ionized electrons are stored."
    )


PICMI_AnyInteraction = PICMI_FieldIonization | PICMI_Interaction
