"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""

from typing import ClassVar, Self, Sequence, get_args, Literal
from pydantic import BaseModel, Field, computed_field, model_validator

from .base import _ClassWithInit, _PICMIModel

class PICMI_BinomialSmoother(_PICMIModel):
    """
    Describes a binomial smoother operator (applied to grids).
    """
    n_pass: Sequence[int] | None = Field(
        default_factory=None,
        description="Vector of integers. Number of passes along each axis"
    )
    compensation: Sequence[bool] | None = Field(
        default=None,
        description="Flags whether to apply compensation along each axis"
    )
    stride: Sequence[int] | None = Field(
        default=None,
        description="Stride along each axis"
    )
    alpha: Sequence[float] | None = Field(
        default=None,
        description="Smoothing coefficients along each axis"
    )

class PICMI_Cartesian1DGrid(_PICMIModel):
    """
    One-dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    Parameters
    ----------
    number_of_cells: vector of integers
        Number of cells along each axis (number of nodes is number_of_cells+1)

    lower_bound: vector of floats
        Position of the node at the lower bound [m]

    upper_bound: vector of floats
        Position of the node at the upper bound [m]

    lower_boundary_conditions: vector of strings
        Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    upper_boundary_conditions: vector of strings
        Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    nx: integer
        Number of cells along X (number of nodes=nx+1)

    xmin: float
        Position of first node along X [m]

    xmax: float
        Position of last node along X [m]

    bc_xmin: vector of strings
        Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_xmax: vector of strings
        Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    moving_window_velocity: vector of floats, optional
        Moving frame velocity [m/s]

    refined_regions: list of lists, optional
        List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
        with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
        lo and hi being vectors of length 2 specifying the extent of the region,
        and refinement_factor defaulting to [2,2] (relative to next lower level)

    lower_bound_particles: vector of floats, optional
        Position of particle lower bound [m]

    upper_bound_particles: vector of floats, optional
        Position of particle upper bound [m]

    xmin_particles: float, optional
        Position of min particle boundary along X [m]

    xmax_particles: float, optional
        Position of max particle boundary along X [m]

    lower_boundary_conditions_particles: vector of strings, optional
        Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal

    upper_boundary_conditions_particles: vector of strings, optional
        Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal

    bc_xmin_particles: string, optional
        Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal

    bc_xmax_particles: string, optional
        Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal

    guard_cells: vector of integers, optional
        Number of guard cells used along each direction

    pml_cells: vector of integers, optional
        Number of Perfectly Matched Layer (PML) cells along each direction

    References
    ----------
    absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
    Based on the Silver-Mueller Radiation Condition, e.g., in

    * A. K. Belhora and L. Pichon, "Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems," 1995,
      https://doi.org/10.1109/20.376322
    * B Engquist and A. Majdat, "Absorbing boundary conditions for numerical simulation of waves," 1977,
      https://doi.org/10.1073/pnas.74.5.1765
    * R. Lehe, "Electromagnetic wave propagation in Particle-In-Cell codes," 2016,
      US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
      https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions: ClassVar[int] = 1

    # Vector forms (the internally-used representation)
    number_of_cells: list[int] | None = None
    lower_bound: list[float] | None = None
    upper_bound: list[float] | None = None
    lower_boundary_conditions: list[str] | None = None
    upper_boundary_conditions: list[str] | None = None
    # Per-axis scalar forms (resolved into the vector forms during validation)
    nx: int | None = None
    xmin: float | None = None
    xmax: float | None = None
    bc_xmin: str | None = None
    bc_xmax: str | None = None
    moving_window_velocity: list[float] | None = None
    refined_regions: list = Field(default_factory=list)
    lower_bound_particles: list[float] | None = None
    upper_bound_particles: list[float] | None = None
    xmin_particles: float | None = None
    xmax_particles: float | None = None
    lower_boundary_conditions_particles: list[str] | None = None
    upper_boundary_conditions_particles: list[str] | None = None
    bc_xmin_particles: str | None = None
    bc_xmax_particles: str | None = None
    guard_cells: list[int] | None = None
    pml_cells: list[int] | None = None

    @model_validator(mode="after")
    def _resolve_grid(self) -> Self:

        # Sanity check and init of input arguments related to grid parameters
        assert (self.number_of_cells is None) and (self.nx is not None) or \
               (self.number_of_cells is not None) and (self.nx is None), \
                'Either number_of_cells or nx must be specified'
        assert (self.lower_bound is None) and (self.xmin is not None) or \
               (self.lower_bound is not None) and (self.xmin is None), \
                'Either lower_bound or xmin must be specified'
        assert (self.upper_bound is None) and (self.xmax is not None) or \
               (self.upper_bound is not None) and (self.xmax is None), \
                'Either upper_bound or xmax must be specified'
        assert (self.lower_boundary_conditions is None) and (self.bc_xmin is not None) or \
               (self.lower_boundary_conditions is not None) and (self.bc_xmin is None), \
                'Either lower_boundary_conditions or bc_xmin'
        assert (self.upper_boundary_conditions is None) and (self.bc_xmax is not None) or \
               (self.upper_boundary_conditions is not None) and (self.bc_xmax is None), \
                'Either upper_boundary_conditions or bc_xmax must be specified'

        if self.number_of_cells is None:
            self.number_of_cells = [self.nx]
        if self.lower_bound is None:
            self.lower_bound = [self.xmin]
        if self.upper_bound is None:
            self.upper_bound = [self.xmax]
        if self.lower_boundary_conditions is None:
            self.lower_boundary_conditions = [self.bc_xmin]
        if self.upper_boundary_conditions is None:
            self.upper_boundary_conditions = [self.bc_xmax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if self.lower_bound_particles is None:
            if (self.xmin_particles is None):
                self.lower_bound_particles = self.lower_bound
            else:
                self.lower_bound_particles = [self.xmin_particles]
        if self.upper_bound_particles is None:
            if (self.xmax_particles is None):
                self.upper_bound_particles = self.upper_bound
            else:
                self.upper_bound_particles = [self.xmax_particles]

        if self.lower_boundary_conditions_particles is None:
            if (self.bc_xmin_particles is None):
                self.lower_boundary_conditions_particles = self.lower_boundary_conditions
            else:
                self.lower_boundary_conditions_particles = [self.bc_xmin_particles]
        if self.upper_boundary_conditions_particles is None:
            if (self.bc_xmax_particles is None):
                self.upper_boundary_conditions_particles = self.upper_boundary_conditions
            else:
                self.upper_boundary_conditions_particles = [self.bc_xmax_particles]

        # Sanity check on dimensionality of vector quantities
        assert len(self.number_of_cells) == 1, 'Wrong number of cells specified'
        assert len(self.lower_bound) == 1, 'Wrong number of lower bounds specified'
        assert len(self.upper_bound) == 1, 'Wrong number of upper bounds specified'
        assert len(self.lower_boundary_conditions) == 1, 'Wrong number of lower boundary conditions specified'
        assert len(self.upper_boundary_conditions) == 1, 'Wrong number of upper boundary conditions specified'
        assert len(self.lower_bound_particles) == 1, 'Wrong number of particle lower bounds specified'
        assert len(self.upper_bound_particles) == 1, 'Wrong number of particle upper bounds specified'
        assert len(self.lower_boundary_conditions_particles) == 1, 'Wrong number of lower particle boundary conditions specified'
        assert len(self.upper_boundary_conditions_particles) == 1, 'Wrong number of upper particle boundary conditions specified'

        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2])
            assert len(region[1]) == 1, 'The lo extent of the refined region must be a vector of length 1'
            assert len(region[2]) == 1, 'The hi extent of the refined region must be a vector of length 1'
            assert len(region[3]) == 1, 'The refinement factor of the refined region must be a vector of length 1'

        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_CylindricalGrid(_PICMIModel):
    """
    Axisymmetric, cylindrical grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    Parameters
    ----------
    number_of_cells: vector of integers
        Number of cells along each axis (number of nodes is number_of_cells+1)

    lower_bound: vector of floats
        Position of the node at the lower bound [m]

    upper_bound: vector of floats
        Position of the node at the upper bound [m]

    lower_boundary_conditions: vector of strings
        Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    upper_boundary_conditions: vector of strings
        Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    nr: integer
        Number of cells along R (number of nodes=nr+1)

    nz: integer
        Number of cells along Z (number of nodes=nz+1)

    n_azimuthal_modes: integer
        Number of azimuthal modes

    rmin: float
        Position of first node along R [m]

    rmax: float
        Position of last node along R [m]

    zmin: float
        Position of first node along Z [m]

    zmax: float
        Position of last node along Z [m]

    bc_rmin: vector of strings
        Boundary condition at min R: One of open, dirichlet, absorbing_silver_mueller, or neumann

    bc_rmax: vector of strings
        Boundary condition at max R: One of open, dirichlet, absorbing_silver_mueller, or neumann

    bc_zmin: vector of strings
        Boundary condition at min Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_zmax: vector of strings
        Boundary condition at max Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    moving_window_velocity: vector of floats, optional
        Moving frame velocity [m/s]

    refined_regions: list of lists, optional
        List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
        with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
        lo and hi being vectors of length 2 specifying the extent of the region,
        and refinement_factor defaulting to [2,2] (relative to next lower level)

    lower_bound_particles: vector of floats, optional
        Position of particle lower bound [m]

    upper_bound_particles: vector of floats, optional
        Position of particle upper bound [m]

    rmin_particles: float, optional
        Position of min particle boundary along R [m]

    rmax_particles: float, optional
        Position of max particle boundary along R [m]

    zmin_particles: float, optional
        Position of min particle boundary along Z [m]

    zmax_particles: float, optional
        Position of max particle boundary along Z [m]

    lower_boundary_conditions_particles: vector of strings, optional
        Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal

    upper_boundary_conditions_particles: vector of strings, optional
        Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal

    bc_rmin_particles: string, optional
        Boundary condition at min R for particles: One of periodic, absorbing, reflect, thermal

    bc_rmax_particles: string, optional
        Boundary condition at max R for particles: One of periodic, absorbing, reflect, thermal

    bc_zmin_particles: string, optional
        Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal

    bc_zmax_particles: string, optional
        Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal

    guard_cells: vector of integers, optional
        Number of guard cells used along each direction

    pml_cells: vector of integers, optional
        Number of Perfectly Matched Layer (PML) cells along each direction

    References
    ----------
    absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
    Based on the Silver-Mueller Radiation Condition, e.g., in

    * A. K. Belhora and L. Pichon, "Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems," 1995,
      https://doi.org/10.1109/20.376322
    * B Engquist and A. Majdat, "Absorbing boundary conditions for numerical simulation of waves," 1977,
      https://doi.org/10.1073/pnas.74.5.1765
    * R. Lehe, "Electromagnetic wave propagation in Particle-In-Cell codes," 2016,
      US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
      https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nr and nz) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions: ClassVar[int] = 2

    # Vector forms (the internally-used representation)
    number_of_cells: list[int] | None = None
    lower_bound: list[float] | None = None
    upper_bound: list[float] | None = None
    lower_boundary_conditions: list[str | None] | None = None
    upper_boundary_conditions: list[str] | None = None
    # Per-axis scalar forms (resolved into the vector forms during validation)
    nr: int | None = None
    nz: int | None = None
    n_azimuthal_modes: int | None = None
    rmin: float | None = None
    rmax: float | None = None
    zmin: float | None = None
    zmax: float | None = None
    bc_rmin: str | None = None
    bc_rmax: str | None = None
    bc_zmin: str | None = None
    bc_zmax: str | None = None
    moving_window_velocity: list[float] | None = None
    refined_regions: list = Field(default_factory=list)
    lower_bound_particles: list[float] | None = None
    upper_bound_particles: list[float] | None = None
    rmin_particles: float | None = None
    rmax_particles: float | None = None
    zmin_particles: float | None = None
    zmax_particles: float | None = None
    lower_boundary_conditions_particles: list[str] | None = None
    upper_boundary_conditions_particles: list[str] | None = None
    bc_rmin_particles: str | None = None
    bc_rmax_particles: str | None = None
    bc_zmin_particles: str | None = None
    bc_zmax_particles: str | None = None
    guard_cells: list[int] | None = None
    pml_cells: list[int] | None = None

    @model_validator(mode="after")
    def _resolve_grid(self) -> Self:

        # Sanity check and init of input arguments related to grid parameters
        assert (self.number_of_cells is None) and (self.nr is not None and self.nz is not None) or \
               (self.number_of_cells is not None) and (self.nr is None and self.nz is None), \
                'Either number_of_cells or nr and nz must be specified'
        assert (self.lower_bound is None) and (self.rmin is not None and self.zmin is not None) or \
               (self.lower_bound is not None) and (self.rmin is None and self.zmin is None), \
                'Either lower_bound or rmin and zmin must be specified'
        assert (self.upper_bound is None) and (self.rmax is not None and self.zmax is not None) or \
               (self.upper_bound is not None) and (self.rmax is None and self.zmax is None), \
                'Either upper_bound or rmax and zmax must be specified'
        # --Allow bc_rmin to be None since it will usually be the axis.
        assert (self.lower_boundary_conditions is None) and (self.bc_zmin is not None) or \
               (self.lower_boundary_conditions is not None) and (self.bc_rmin is None and self.bc_zmin is None), \
                'Either lower_boundary_conditions or bc_rmin and bc_zmin must be specified'
        assert (self.upper_boundary_conditions is None) and (self.bc_rmax is not None and self.bc_zmax is not None) or \
               (self.upper_boundary_conditions is not None) and (self.bc_rmax is None and self.bc_zmax is None), \
                'Either upper_boundary_conditions or bc_rmax and bc_zmax must be specified'

        if self.number_of_cells is None:
            self.number_of_cells = [self.nr, self.nz]
        if self.lower_bound is None:
            self.lower_bound = [self.rmin, self.zmin]
        if self.upper_bound is None:
            self.upper_bound = [self.rmax, self.zmax]
        if self.lower_boundary_conditions is None:
            self.lower_boundary_conditions = [self.bc_rmin, self.bc_zmin]
        if self.upper_boundary_conditions is None:
            self.upper_boundary_conditions = [self.bc_rmax, self.bc_zmax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if self.lower_bound_particles is None:
            if (self.rmin_particles is None) and (self.zmin_particles is None):
                self.lower_bound_particles = self.lower_bound
            else:
                self.lower_bound_particles = [self.rmin_particles, self.zmin_particles]
        if self.upper_bound_particles is None:
            if (self.rmax_particles is None) and (self.zmax_particles is None):
                self.upper_bound_particles = self.upper_bound
            else:
                self.upper_bound_particles = [self.rmax_particles, self.zmax_particles]

        if self.lower_boundary_conditions_particles is None:
            if (self.bc_rmin_particles is None) and (self.bc_zmin_particles is None):
                self.lower_boundary_conditions_particles = self.lower_boundary_conditions
            else:
                self.lower_boundary_conditions_particles = [self.bc_rmin_particles, self.bc_zmin_particles]
        if self.upper_boundary_conditions_particles is None:
            if (self.bc_rmax_particles is None) and (self.bc_zmax_particles is None):
                self.upper_boundary_conditions_particles = self.upper_boundary_conditions
            else:
                self.upper_boundary_conditions_particles = [self.bc_rmax_particles, self.bc_zmax_particles]

        # Sanity check on dimensionality of vector quantities
        assert len(self.number_of_cells) == 2, 'Wrong number of cells specified'
        assert len(self.lower_bound) == 2, 'Wrong number of lower bounds specified'
        assert len(self.upper_bound) == 2, 'Wrong number of upper bounds specified'
        assert len(self.lower_boundary_conditions) == 2, 'Wrong number of lower boundary conditions specified'
        assert len(self.upper_boundary_conditions) == 2, 'Wrong number of upper boundary conditions specified'

        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2,2])
            assert len(region[1]) == 2, 'The lo extent of the refined region must be a vector of length 2'
            assert len(region[2]) == 2, 'The hi extent of the refined region must be a vector of length 2'
            assert len(region[3]) == 2, 'The refinement factor of the refined region must be a vector of length 2'

        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2,2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_Cartesian2DGrid(_PICMIModel):
    """
    Two dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    Parameters
    ----------
    number_of_cells: vector of integers
        Number of cells along each axis (number of nodes is number_of_cells+1)

    lower_bound: vector of floats
        Position of the node at the lower bound [m]

    upper_bound: vector of floats
        Position of the node at the upper bound [m]

    lower_boundary_conditions: vector of strings
        Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    upper_boundary_conditions: vector of strings
        Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    nx: integer
        Number of cells along X (number of nodes=nx+1)

    ny: integer
        Number of cells along Y (number of nodes=ny+1)

    xmin: float
        Position of first node along X [m]

    xmax: float
        Position of last node along X [m]

    ymin: float
        Position of first node along Y [m]

    ymax: float
        Position of last node along Y [m]

    bc_xmin: vector of strings
        Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_xmax: vector of strings
        Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_ymin: vector of strings
        Boundary condition at min Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_ymax: vector of strings
        Boundary condition at max Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    moving_window_velocity: vector of floats, optional
        Moving frame velocity [m/s]

    refined_regions: list of lists, optional
        List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
        with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
        lo and hi being vectors of length 2 specifying the extent of the region,
        and refinement_factor defaulting to [2,2] (relative to next lower level)

    lower_bound_particles: vector of floats, optional
        Position of particle lower bound [m]

    upper_bound_particles: vector of floats, optional
        Position of particle upper bound [m]

    xmin_particles: float, optional
        Position of min particle boundary along X [m]

    xmax_particles: float, optional
        Position of max particle boundary along X [m]

    ymin_particles: float, optional
        Position of min particle boundary along Y [m]

    ymax_particles: float, optional
        Position of max particle boundary along Y [m]

    lower_boundary_conditions_particles: vector of strings, optional
        Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal

    upper_boundary_conditions_particles: vector of strings, optional
        Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal

    bc_xmin_particles: string, optional
        Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal

    bc_xmax_particles: string, optional
        Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal

    bc_ymin_particles: string, optional
        Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal

    bc_ymax_particles: string, optional
        Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal

    guard_cells: vector of integers, optional
        Number of guard cells used along each direction

    pml_cells: vector of integers, optional
        Number of Perfectly Matched Layer (PML) cells along each direction

    References
    ----------
    absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
    Based on the Silver-Mueller Radiation Condition, e.g., in

    * A. K. Belhora and L. Pichon, "Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems," 1995,
      https://doi.org/10.1109/20.376322
    * B Engquist and A. Majdat, "Absorbing boundary conditions for numerical simulation of waves," 1977,
      https://doi.org/10.1073/pnas.74.5.1765
    * R. Lehe, "Electromagnetic wave propagation in Particle-In-Cell codes," 2016,
      US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
      https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx and ny) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions: ClassVar[int] = 2

    # Vector forms (the internally-used representation)
    number_of_cells: list[int] | None = None
    lower_bound: list[float] | None = None
    upper_bound: list[float] | None = None
    lower_boundary_conditions: list[str] | None = None
    upper_boundary_conditions: list[str] | None = None
    # Per-axis scalar forms (resolved into the vector forms during validation)
    nx: int | None = None
    ny: int | None = None
    xmin: float | None = None
    xmax: float | None = None
    ymin: float | None = None
    ymax: float | None = None
    bc_xmin: str | None = None
    bc_xmax: str | None = None
    bc_ymin: str | None = None
    bc_ymax: str | None = None
    moving_window_velocity: list[float] | None = None
    refined_regions: list = Field(default_factory=list)
    lower_bound_particles: list[float] | None = None
    upper_bound_particles: list[float] | None = None
    xmin_particles: float | None = None
    xmax_particles: float | None = None
    ymin_particles: float | None = None
    ymax_particles: float | None = None
    lower_boundary_conditions_particles: list[str] | None = None
    upper_boundary_conditions_particles: list[str] | None = None
    bc_xmin_particles: str | None = None
    bc_xmax_particles: str | None = None
    bc_ymin_particles: str | None = None
    bc_ymax_particles: str | None = None
    guard_cells: list[int] | None = None
    pml_cells: list[int] | None = None

    @model_validator(mode="after")
    def _resolve_grid(self) -> Self:

        # Sanity check and init of input arguments related to grid parameters
        assert (self.number_of_cells is None) and (self.nx is not None and self.ny is not None) or \
               (self.number_of_cells is not None) and (self.nx is None and self.ny is None), \
                'Either number_of_cells or nx and ny must be specified'
        assert (self.lower_bound is None) and (self.xmin is not None and self.ymin is not None) or \
               (self.lower_bound is not None) and (self.xmin is None and self.ymin is None), \
                'Either lower_bound or xmin and ymin must be specified'
        assert (self.upper_bound is None) and (self.xmax is not None and self.ymax is not None) or \
               (self.upper_bound is not None) and (self.xmax is None and self.ymax is None), \
                'Either upper_bound or xmax and ymax must be specified'
        assert (self.lower_boundary_conditions is None) and (self.bc_xmin is not None and self.bc_ymin is not None) or \
               (self.lower_boundary_conditions is not None) and (self.bc_xmin is None and self.bc_ymin is None), \
                'Either lower_boundary_conditions or bc_xmin and bc_ymin must be specified'
        assert (self.upper_boundary_conditions is None) and (self.bc_xmax is not None and self.bc_ymax is not None) or \
               (self.upper_boundary_conditions is not None) and (self.bc_xmax is None and self.bc_ymax is None), \
                'Either upper_boundary_conditions or bc_xmax and bc_ymax must be specified'

        if self.number_of_cells is None:
            self.number_of_cells = [self.nx, self.ny]
        if self.lower_bound is None:
            self.lower_bound = [self.xmin, self.ymin]
        if self.upper_bound is None:
            self.upper_bound = [self.xmax, self.ymax]
        if self.lower_boundary_conditions is None:
            self.lower_boundary_conditions = [self.bc_xmin, self.bc_ymin]
        if self.upper_boundary_conditions is None:
            self.upper_boundary_conditions = [self.bc_xmax, self.bc_ymax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if self.lower_bound_particles is None:
            if (self.xmin_particles is None) and (self.ymin_particles is None):
                self.lower_bound_particles = self.lower_bound
            else:
                self.lower_bound_particles = [self.xmin_particles, self.ymin_particles]
        if self.upper_bound_particles is None:
            if (self.xmax_particles is None) and (self.ymax_particles is None):
                self.upper_bound_particles = self.upper_bound
            else:
                self.upper_bound_particles = [self.xmax_particles, self.ymax_particles]

        if self.lower_boundary_conditions_particles is None:
            if (self.bc_xmin_particles is None) and (self.bc_ymin_particles is None):
                self.lower_boundary_conditions_particles = self.lower_boundary_conditions
            else:
                self.lower_boundary_conditions_particles = [self.bc_xmin_particles, self.bc_ymin_particles]
        if self.upper_boundary_conditions_particles is None:
            if (self.bc_xmax_particles is None) and (self.bc_ymax_particles is None):
                self.upper_boundary_conditions_particles = self.upper_boundary_conditions
            else:
                self.upper_boundary_conditions_particles = [self.bc_xmax_particles, self.bc_ymax_particles]

        # Sanity check on dimensionality of vector quantities
        assert len(self.number_of_cells) == 2, 'Wrong number of cells specified'
        assert len(self.lower_bound) == 2, 'Wrong number of lower bounds specified'
        assert len(self.upper_bound) == 2, 'Wrong number of upper bounds specified'
        assert len(self.lower_boundary_conditions) == 2, 'Wrong number of lower boundary conditions specified'
        assert len(self.upper_boundary_conditions) == 2, 'Wrong number of upper boundary conditions specified'
        assert len(self.lower_bound_particles) == 2, 'Wrong number of particle lower bounds specified'
        assert len(self.upper_bound_particles) == 2, 'Wrong number of particle upper bounds specified'
        assert len(self.lower_boundary_conditions_particles) == 2, 'Wrong number of lower particle boundary conditions specified'
        assert len(self.upper_boundary_conditions_particles) == 2, 'Wrong number of upper particle boundary conditions specified'

        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2,2])
            assert len(region[1]) == 2, 'The lo extent of the refined region must be a vector of length 2'
            assert len(region[2]) == 2, 'The hi extent of the refined region must be a vector of length 2'
            assert len(region[3]) == 2, 'The refinement factor of the refined region must be a vector of length 2'

        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2,2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_Cartesian3DGrid(_PICMIModel):
    """
    Three dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    Parameters
    ----------
    number_of_cells: vector of integers
        Number of cells along each axis (number of nodes is number_of_cells+1)

    lower_bound: vector of floats
        Position of the node at the lower bound [m]

    upper_bound: vector of floats
        Position of the node at the upper bound [m]

    lower_boundary_conditions: vector of strings
        Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    upper_boundary_conditions: vector of strings
        Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    nx: integer
        Number of cells along X (number of nodes=nx+1)

    ny: integer
        Number of cells along Y (number of nodes=ny+1)

    nz: integer
        Number of cells along Z (number of nodes=nz+1)

    xmin: float
        Position of first node along X [m]

    xmax: float
        Position of last node along X [m]

    ymin: float
        Position of first node along Y [m]

    ymax: float
        Position of last node along Y [m]

    zmin: float
        Position of first node along Z [m]

    zmax: float
        Position of last node along Z [m]

    bc_xmin: string
        Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_xmax: string
        Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_ymin: string
        Boundary condition at min Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_ymax: string
        Boundary condition at max Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_zmin: string
        Boundary condition at min Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    bc_zmax: string
        Boundary condition at max Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann

    moving_window_velocity: vector of floats, optional
        Moving frame velocity [m/s]

    refined_regions: list of lists, optional
        List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
        with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
        lo and hi being vectors of length 3 specifying the extent of the region,
        and refinement_factor defaulting to [2,2,2] (relative to next lower level)

    lower_bound_particles: vector of floats, optional
        Position of particle lower bound [m]

    upper_bound_particles: vector of floats, optional
        Position of particle upper bound [m]

    xmin_particles: float, optional
        Position of min particle boundary along X [m]

    xmax_particles: float, optional
        Position of max particle boundary along X [m]

    ymin_particles: float, optional
        Position of min particle boundary along Y [m]

    ymax_particles: float, optional
        Position of max particle boundary along Y [m]

    zmin_particles float, optional
        Position of min particle boundary along Z [m]

    zmax_particles: float, optional
        Position of max particle boundary along Z [m]

    lower_boundary_conditions_particles: vector of strings, optional
        Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal

    upper_boundary_conditions_particles: vector of strings, optional
        Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal

    bc_xmin_particles: string, optional
        Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal

    bc_xmax_particles: string, optional
        Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal

    bc_ymin_particles: string, optional
        Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal

    bc_ymax_particles: string, optional
        Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal

    bc_zmin_particles: string, optional
        Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal

    bc_zmax_particles: string, optional
        Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal

    guard_cells: vector of integers, optional
        Number of guard cells used along each direction

    pml_cells: vector of integers, optional
        Number of Perfectly Matched Layer (PML) cells along each direction

    References
    ----------
    absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
    Based on the Silver-Mueller Radiation Condition, e.g., in

    * A. K. Belhora and L. Pichon, "Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems," 1995,
      https://doi.org/10.1109/20.376322
    * B Engquist and A. Majdat, "Absorbing boundary conditions for numerical simulation of waves," 1977,
      https://doi.org/10.1073/pnas.74.5.1765
    * R. Lehe, "Electromagnetic wave propagation in Particle-In-Cell codes," 2016,
      US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
      https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx, ny, and nz) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions: ClassVar[int] = 3

    # Vector forms (the internally-used representation)
    number_of_cells: list[int] | None = None
    lower_bound: list[float] | None = None
    upper_bound: list[float] | None = None
    lower_boundary_conditions: list[str] | None = None
    upper_boundary_conditions: list[str] | None = None
    # Per-axis scalar forms (resolved into the vector forms during validation)
    nx: int | None = None
    ny: int | None = None
    nz: int | None = None
    xmin: float | None = None
    xmax: float | None = None
    ymin: float | None = None
    ymax: float | None = None
    zmin: float | None = None
    zmax: float | None = None
    bc_xmin: str | None = None
    bc_xmax: str | None = None
    bc_ymin: str | None = None
    bc_ymax: str | None = None
    bc_zmin: str | None = None
    bc_zmax: str | None = None
    moving_window_velocity: list[float] | None = None
    refined_regions: list = Field(default_factory=list)
    lower_bound_particles: list[float] | None = None
    upper_bound_particles: list[float] | None = None
    xmin_particles: float | None = None
    xmax_particles: float | None = None
    ymin_particles: float | None = None
    ymax_particles: float | None = None
    zmin_particles: float | None = None
    zmax_particles: float | None = None
    lower_boundary_conditions_particles: list[str] | None = None
    upper_boundary_conditions_particles: list[str] | None = None
    bc_xmin_particles: str | None = None
    bc_xmax_particles: str | None = None
    bc_ymin_particles: str | None = None
    bc_ymax_particles: str | None = None
    bc_zmin_particles: str | None = None
    bc_zmax_particles: str | None = None
    guard_cells: list[int] | None = None
    pml_cells: list[int] | None = None

    @model_validator(mode="after")
    def _resolve_grid(self) -> Self:

        # Sanity check and init of input arguments related to grid parameters
        assert (self.number_of_cells is None) and (self.nx is not None and self.ny is not None and self.nz is not None) or \
               (self.number_of_cells is not None) and (self.nx is None and self.ny is None and self.nz is None), \
                'Either number_of_cells or nx, ny, and nz must be specified'
        assert (self.lower_bound is None) and (self.xmin is not None and self.ymin is not None and self.zmin is not None) or \
               (self.lower_bound is not None) and (self.xmin is None and self.ymin is None and self.zmin is None), \
                'Either lower_bound or xmin, ymin, and zmin must be specified'
        assert (self.upper_bound is None) and (self.xmax is not None and self.ymax is not None and self.zmax is not None) or \
               (self.upper_bound is not None) and (self.xmax is None and self.ymax is None and self.zmax is None), \
                'Either upper_bound or xmax, ymax, and zmax must be specified'
        assert (self.lower_boundary_conditions is None) and (self.bc_xmin is not None and self.bc_ymin is not None and self.bc_zmin is not None) or \
               (self.lower_boundary_conditions is not None) and (self.bc_xmin is None and self.bc_ymin is None and self.bc_zmin is None), \
                'Either lower_boundary_conditions or bc_xmin, bc_ymin, and bc_zmin must be specified'
        assert (self.upper_boundary_conditions is None) and (self.bc_xmax is not None and self.bc_ymax is not None and self.bc_zmax is not None) or \
               (self.upper_boundary_conditions is not None) and (self.bc_xmax is None and self.bc_ymax is None and self.bc_zmax is None), \
                'Either upper_boundary_conditions or bc_xmax, bc_ymax, and bc_zmax must be specified'

        if self.number_of_cells is None:
            self.number_of_cells = [self.nx, self.ny, self.nz]
        if self.lower_bound is None:
            self.lower_bound = [self.xmin, self.ymin, self.zmin]
        if self.upper_bound is None:
            self.upper_bound = [self.xmax, self.ymax, self.zmax]
        if self.lower_boundary_conditions is None:
            self.lower_boundary_conditions = [self.bc_xmin, self.bc_ymin, self.bc_zmin]
        if self.upper_boundary_conditions is None:
            self.upper_boundary_conditions = [self.bc_xmax, self.bc_ymax, self.bc_zmax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if self.lower_bound_particles is None:
            if (self.xmin_particles is None) and (self.ymin_particles is None) and (self.zmin_particles is None):
                self.lower_bound_particles = self.lower_bound
            else:
                self.lower_bound_particles = [self.xmin_particles, self.ymin_particles, self.zmin_particles]
        if self.upper_bound_particles is None:
            if (self.xmax_particles is None) and (self.ymax_particles is None) and (self.zmax_particles is None):
                self.upper_bound_particles = self.upper_bound
            else:
                self.upper_bound_particles = [self.xmax_particles, self.ymax_particles, self.zmax_particles]

        if self.lower_boundary_conditions_particles is None:
            if (self.bc_xmin_particles is None) and (self.bc_ymin_particles is None) and (self.bc_zmin_particles is None):
                self.lower_boundary_conditions_particles = self.lower_boundary_conditions
            else:
                self.lower_boundary_conditions_particles = [self.bc_xmin_particles, self.bc_ymin_particles, self.bc_zmin_particles]
        if self.upper_boundary_conditions_particles is None:
            if (self.bc_xmax_particles is None) and (self.bc_ymax_particles is None) and (self.bc_zmax_particles is None):
                self.upper_boundary_conditions_particles = self.upper_boundary_conditions
            else:
                self.upper_boundary_conditions_particles = [self.bc_xmax_particles, self.bc_ymax_particles, self.bc_zmax_particles]

        # Sanity check on number of arguments of vector quantities
        assert len(self.number_of_cells) == 3, 'Wrong number of cells specified'
        assert len(self.lower_bound) == 3, 'Wrong number of lower bounds specified'
        assert len(self.upper_bound) == 3, 'Wrong number of upper bounds specified'
        assert len(self.lower_boundary_conditions) == 3, 'Wrong number of lower boundary conditions specified'
        assert len(self.upper_boundary_conditions) == 3, 'Wrong number of upper boundary conditions specified'
        assert len(self.lower_bound_particles) == 3, 'Wrong number of particle lower bounds specified'
        assert len(self.upper_bound_particles) == 3, 'Wrong number of particle upper bounds specified'
        assert len(self.lower_boundary_conditions_particles) == 3, 'Wrong number of particle lower boundary conditions specified'
        assert len(self.upper_boundary_conditions_particles) == 3, 'Wrong number of particle upper boundary conditions specified'

        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2,2,2])
            assert len(region[1]) == 3, 'The lo extent of the refined region must be a vector of length 3'
            assert len(region[2]) == 3, 'The hi extent of the refined region must be a vector of length 3'
            assert len(region[3]) == 3, 'The refinement factor of the refined region must be a vector of length 3'

        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2,2,2]):
        """Add a refined region.

        Parameters
        ----------
        level: integer
            The refinement level, with 1 being the first level of refinement, 2 being the second etc.

        lo, hi: vectors of floats
            Each is a vector of length 3 specifying the extent of the region

        refinement_factor: vector of integers, optional
            Defaulting to [2,2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


PICMI_AnyGrid = PICMI_CylindricalGrid | PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid

class PICMI_ElectromagneticSolver(_PICMIModel):
    """
    Electromagnetic field solver.

    The advance method used to solve Maxwell's equations. The default method is code dependent.

    Method options:

    - 'Yee': standard solver using the staggered Yee grid (https://doi.org/10.1109/TAP.1966.1138693)
    - 'CKC': solver with the extended Cole-Karkkainen-Cowan stencil with better dispersion properties (https://doi.org/10.1103/PhysRevSTAB.16.041303)
    - 'Lehe': CKC-style solver with modified dispersion (https://doi.org/10.1103/PhysRevSTAB.16.021301)
    - 'PSTD': Spectral solver with finite difference in time domain, e.g., Q. H. Liu, Letters 15 (3) (1997) 158–165
    - 'PSATD': Spectral solver with analytic in time domain (https://doi.org/10.1016/j.jcp.2013.03.010)
    - 'DS': Directional Splitting after Yasuhiko Sentoku (https://doi.org/10.1140/epjd/e2014-50162-y)
    - 'ECT': Enlarged Cell Technique solver, allowing internal conductors (https://doi.org/10.1109/APS.2005.1551259)
    """

    @computed_field
    def methods_list(self) -> list[str]:
        # Retained for backwards compatibility reasons.
        # The type annotation of `method` is the ground-truth.
        return list(get_args(type(self).__annotations__['method']))

    grid: PICMI_AnyGrid = Field(description="Grid object for the diagnostic")
    method: Literal['Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD', 'DS', 'ECT'] | None = Field(
        default=None,
        description="The advance method use to solve Maxwell's equations. The default method is code dependent."
    )
    stencil_order: Sequence[int] | None = Field(
        default=None,
        description="Order of stencil for each axis (-1=infinite)"
    )
    cfl: float | None = Field(
        default=None,
        description="Fraction of the Courant-Friedrich-Lewy criteria [1]"
    )
    source_smoother: PICMI_BinomialSmoother | None = Field(
        default=None,
        description="Smoother object to apply to the sources"
    )
    field_smoother: PICMI_BinomialSmoother | None = Field(
        default=None,
        description="Smoother object to apply to the fields"
    )
    subcycling: int | None = Field(
        default=None,
        description="Level of subcycling for the GPSTD solver"
    )
    galilean_velocity: Sequence[float] | None = Field(
        default=None,
        description="Velocity of Galilean reference frame [m/s]"
    )
    divE_cleaning: bool | None = Field(
        default=None,
        description="Solver uses div(E) cleaning if True"
    )
    divB_cleaning: bool | None = Field(
        default=None,
        description="Solver uses div(B) cleaning if True"
    )
    pml_divE_cleaning: bool | None = Field(
        default=None,
        description="Solver uses div(E) cleaning in the PML if True"
    )
    pml_divB_cleaning: bool | None = Field(
        default=None,
        description="Solver uses div(B) cleaning in the PML if True"
    )


class PICMI_ElectrostaticSolver(_ClassWithInit):
    """
    Electrostatic field solver

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    method: string
        One of 'FFT', or 'Multigrid'

    required_precision: float, optional
        Level of precision required for iterative solvers

    maximum_iterations: integer, optional
        Maximum number of iterations for iterative solvers
    """

    methods_list = ['FFT', 'Multigrid']

    def __init__(self, grid, method=None,
                 required_precision=None, maximum_iterations=None, **kw):

        assert method is None or method in PICMI_ElectrostaticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_ElectrostaticSolver.methods_list))

        self.grid = grid
        self.method = method
        self.required_precision = required_precision
        self.maximum_iterations = maximum_iterations

        self.handle_init(kw)


class PICMI_MagnetostaticSolver(_ClassWithInit):
    """
    Magnetostatic field solver

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    method: string
        One of 'FFT', or 'Multigrid'
    """

    methods_list = ['FFT', 'Multigrid']

    def __init__(self, grid, method=None, **kw):

        assert method is None or method in PICMI_MagnetostaticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_MagnetostaticSolver.methods_list))

        self.grid = grid
        self.method = method

        self.handle_init(kw)
