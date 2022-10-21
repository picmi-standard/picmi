"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in this file are related to interactions (e.g. field ionization, collisions, QED)
"""

from .base import _ClassWithInit

class PICMI_FieldIonization(_ClassWithInit):
    """
    Field ionization on an ion species

    Parameters
    ----------
    model: string
        Ionization model, e.g. "ADK"

    ionized_species: species instance
        Species that is ionized

    product_species: species instance
        Species in which ionized electrons are stored.
    """
    def __init__(self, model, ionized_species, product_species, **kw):
        self.model = model
        self.ionized_species = ionized_species
        self.product_species = product_species

        self.handle_init(kw)

