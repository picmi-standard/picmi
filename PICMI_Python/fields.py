"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""

from typing import ClassVar, Literal, NamedTuple, Self, get_args
from pydantic import Field, PrivateAttr, model_validator

from .base import _PICMIModel, PICMI_Grid, PICMI_Solver, resolve_once


class _AxisGroup(NamedTuple):
    """A vector parameter of a grid and its per-axis forms, e.g., number_of_cells and nx, ny, nz"""

    vector: str
    per_axis: tuple[str, ...]
    # the vector that the axes, which are not specified, default to
    default: str | None = None
    # whether None is a valid value of an axis (e.g., the boundary condition on the axis in RZ)
    none_is_value: bool = False


class _AxisGroupState(NamedTuple):
    values: list
    per_axis: list
    defaulted_axes: frozenset


class _PICMIGrid(PICMI_Grid):
    # Base of the grids (without docstring, so that it is not prepended to the grids' ones).
    #
    # Grid parameters can be specified as vectors (e.g., number_of_cells) or per axis (e.g.,
    # nx, ny, nz). Both forms are kept in sync: at construction, the vector is used if both are
    # given. Afterwards, assigning to either form updates the other one, and assigning None
    # restores the current value (or for the particle boundaries, returns to the default).
    # Particle boundaries that are not specified default to the field boundaries, and follow
    # them when those change later on.

    # The vector parameters with their per-axis forms, defined by each grid. Groups with a
    # default are listed after the group of their default vector.
    _axis_groups: ClassVar[tuple[_AxisGroup, ...]] = ()

    # Per vector parameter: the state as of the last validation (None before the first one).
    _axis_state: dict[str, _AxisGroupState] | None = PrivateAttr(default=None)

    def _resolve_axis_groups(self):
        """Resolve and synchronize the vector and per-axis forms of all grid parameters.

        Validation re-runs on every assignment and when the grid is passed to another PICMI
        object, so this compares against the state of the last validation to find out which
        form was assigned.
        """
        new_state = {}
        for group in self._axis_groups:
            vector = getattr(self, group.vector)
            per_axis = [getattr(self, name) for name in group.per_axis]
            previous = None if self._axis_state is None else self._axis_state[group.vector]

            if previous is None:
                # initial input: the vector is used if given
                if vector is not None:
                    values, defaulted_axes = list(vector), set()
                else:
                    values = per_axis
                    defaulted_axes = set()
                    if group.default is not None:
                        defaulted_axes = {axis for axis, value in enumerate(values) if value is None}
            elif vector != previous.values:
                if vector is not None:
                    values, defaulted_axes = list(vector), set()
                elif group.default is not None:
                    values, defaulted_axes = list(previous.values), set(range(len(previous.values)))
                else:
                    values, defaulted_axes = list(previous.values), set(previous.defaulted_axes)
            else:
                values, defaulted_axes = list(previous.values), set(previous.defaulted_axes)
                for axis, (value, previous_value) in enumerate(zip(per_axis, previous.per_axis)):
                    if value == previous_value:
                        continue
                    if value is not None or group.none_is_value:
                        values[axis] = value
                        defaulted_axes.discard(axis)
                    elif group.default is not None:
                        defaulted_axes.add(axis)

            if defaulted_axes:
                default = getattr(self, group.default)
                for axis in defaulted_axes:
                    if axis < len(values) and default is not None and axis < len(default):
                        values[axis] = default[axis]

            if vector != values:
                setattr(self, group.vector, values)
            # a vector of the wrong length is reported by the dimensionality checks of the grid
            if len(values) == len(group.per_axis):
                for name, value in zip(group.per_axis, values):
                    if getattr(self, name) != value:
                        setattr(self, name, value)
            new_state[group.vector] = _AxisGroupState(
                list(getattr(self, group.vector)),
                [getattr(self, name) for name in group.per_axis],
                frozenset(defaulted_axes),
            )
        self._axis_state = new_state

    def _resolve_cartesian_grid(self):
        """Check and resolve the parameters of a Cartesian grid

        The Cartesian grids differ only in their axes, so they share these checks. The
        cylindrical grid checks its parameters itself, since its radial axis is special.
        """
        for group in self._axis_groups:
            if group.default is not None or getattr(self, group.vector) is not None:
                # given as a vector, or it follows another parameter if it is not given
                continue
            axes = group.per_axis
            listed = (
                " and ".join(axes)
                if len(axes) < 3
                else f"{', '.join(axes[:-1])}, and {axes[-1]}"
            )
            assert all(getattr(self, axis) is not None for axis in axes), (
                f"Either {group.vector} or {listed} must be specified"
            )

        # Resolve and synchronize the vector and per-axis forms, see above.
        # By default, if not specified, the particle boundary values and conditions are the
        # same as the field boundary values and conditions.
        self._resolve_axis_groups()

        # Sanity check on dimensionality of vector quantities
        for group in self._axis_groups:
            values = getattr(self, group.vector)
            assert len(values) == self.number_of_dimensions, (
                f"Wrong number of values in {group.vector}: {len(values)} "
                f"instead of {self.number_of_dimensions}"
            )

        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2] * self.number_of_dimensions)
            for index, name in (
                (1, "lo extent"),
                (2, "hi extent"),
                (3, "refinement factor"),
            ):
                assert len(region[index]) == self.number_of_dimensions, (
                    f"The {name} of the refined region must be a vector of length "
                    f"{self.number_of_dimensions}"
                )


def _grid_axis_groups(axes, lower_boundary_condition_can_be_none=False):
    """The axis groups of a grid with the given axis names, e.g., ("x", "y", "z")"""
    return (
        _AxisGroup("number_of_cells", tuple(f"n{axis}" for axis in axes)),
        _AxisGroup("lower_bound", tuple(f"{axis}min" for axis in axes)),
        _AxisGroup("upper_bound", tuple(f"{axis}max" for axis in axes)),
        _AxisGroup(
            "lower_boundary_conditions",
            tuple(f"bc_{axis}min" for axis in axes),
            none_is_value=lower_boundary_condition_can_be_none,
        ),
        _AxisGroup("upper_boundary_conditions", tuple(f"bc_{axis}max" for axis in axes)),
        _AxisGroup(
            "lower_bound_particles",
            tuple(f"{axis}min_particles" for axis in axes),
            default="lower_bound",
        ),
        _AxisGroup(
            "upper_bound_particles",
            tuple(f"{axis}max_particles" for axis in axes),
            default="upper_bound",
        ),
        _AxisGroup(
            "lower_boundary_conditions_particles",
            tuple(f"bc_{axis}min_particles" for axis in axes),
            default="lower_boundary_conditions",
        ),
        _AxisGroup(
            "upper_boundary_conditions_particles",
            tuple(f"bc_{axis}max_particles" for axis in axes),
            default="upper_boundary_conditions",
        ),
    )


class PICMI_BinomialSmoother(_PICMIModel):
    """
    Describes a binomial smoother operator (applied to grids).
    """

    n_pass: list[int] | None = Field(
        default=None,
        description="Vector of integers. Number of passes along each axis",
    )
    compensation: list[bool] | None = Field(
        default=None, description="Flags whether to apply compensation along each axis"
    )
    stride: list[int] | None = Field(
        default=None, description="Stride along each axis"
    )
    alpha: list[float] | None = Field(
        default=None, description="Smoothing coefficients along each axis"
    )


class PICMI_Cartesian1DGrid(_PICMIGrid):
    """
    One-dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

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
    # vector of values (such as number_of_cells). Both forms are kept in sync (see _PICMIGrid), but the
    # implementation should use the vectors to access the user input.

    number_of_dimensions: ClassVar[int] = 1
    _axis_groups: ClassVar[tuple[_AxisGroup, ...]] = _grid_axis_groups(("x",))

    # Vector forms (the internally-used representation)
    number_of_cells: list[int] | None = Field(
        default=None,
        description="Number of cells along each axis (number of nodes is number_of_cells+1)",
    )
    lower_bound: list[float] | None = Field(
        default=None, description="Position of the node at the lower bound [m]"
    )
    upper_bound: list[float] | None = Field(
        default=None, description="Position of the node at the upper bound [m]"
    )
    lower_boundary_conditions: list[str | None] | None = Field(
        default=None,
        description="Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    upper_boundary_conditions: list[str | None] | None = Field(
        default=None,
        description="Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    # Per-axis scalar forms (resolved into the vector forms during validation)
    nx: int | None = Field(
        default=None, description="Number of cells along X (number of nodes=nx+1)"
    )
    xmin: float | None = Field(
        default=None, description="Position of first node along X [m]"
    )
    xmax: float | None = Field(
        default=None, description="Position of last node along X [m]"
    )
    bc_xmin: str | None = Field(
        default=None,
        description="Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_xmax: str | None = Field(
        default=None,
        description="Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    moving_window_velocity: list[float] | None = Field(
        default=None, description="Moving frame velocity [m/s]"
    )
    refined_regions: list = Field(
        default_factory=list,
        description="List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor], with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc, lo and hi being vectors of length 1 specifying the extent of the region, and refinement_factor defaulting to [2] (relative to next lower level)",
    )
    lower_bound_particles: list[float] | None = Field(
        default=None, description="Position of particle lower bound [m]"
    )
    upper_bound_particles: list[float] | None = Field(
        default=None, description="Position of particle upper bound [m]"
    )
    xmin_particles: float | None = Field(
        default=None, description="Position of min particle boundary along X [m]"
    )
    xmax_particles: float | None = Field(
        default=None, description="Position of max particle boundary along X [m]"
    )
    lower_boundary_conditions_particles: list[str | None] | None = Field(
        default=None,
        description="Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal",
    )
    upper_boundary_conditions_particles: list[str | None] | None = Field(
        default=None,
        description="Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal",
    )
    bc_xmin_particles: str | None = Field(
        default=None,
        description="Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_xmax_particles: str | None = Field(
        default=None,
        description="Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal",
    )
    guard_cells: list[int] | None = Field(
        default=None, description="Number of guard cells used along each direction"
    )
    pml_cells: list[int] | None = Field(
        default=None,
        description="Number of Perfectly Matched Layer (PML) cells along each direction",
    )

    @model_validator(mode="after")
    @resolve_once
    def _resolve_grid(self) -> Self:
        self._resolve_cartesian_grid()
        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2]):
        """Add a refined region.

        Parameters
        ----------
        level : int
            The refinement level, with 1 being the first level of refinement, 2 being the second etc.

        lo, hi : list of float
            The lower (lo) and upper (hi) corner of the region, each a vector of length 1

        refinement_factor : list of int, optional
            The refinement factor relative to the next lower level (default: [2])
        """
        # assign (instead of appending in place), so that the region is validated and marked as set
        self.refined_regions = [*self.refined_regions, [level, lo, hi, refinement_factor]]


class PICMI_CylindricalGrid(_PICMIGrid):
    """
    Axisymmetric, cylindrical grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

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
    # vector of values (such as number_of_cells). Both forms are kept in sync (see _PICMIGrid), but the
    # implementation should use the vectors to access the user input.

    number_of_dimensions: ClassVar[int] = 2
    _axis_groups: ClassVar[tuple[_AxisGroup, ...]] = _grid_axis_groups(("r", "z"), lower_boundary_condition_can_be_none=True)

    # Vector forms (the internally-used representation)
    number_of_cells: list[int] | None = Field(
        default=None,
        description="Number of cells along each axis (number of nodes is number_of_cells+1)",
    )
    lower_bound: list[float] | None = Field(
        default=None, description="Position of the node at the lower bound [m]"
    )
    upper_bound: list[float] | None = Field(
        default=None, description="Position of the node at the upper bound [m]"
    )
    lower_boundary_conditions: list[str | None] | None = Field(
        default=None,
        description="Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    upper_boundary_conditions: list[str | None] | None = Field(
        default=None,
        description="Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    # Per-axis scalar forms (resolved into the vector forms during validation)
    nr: int | None = Field(
        default=None, description="Number of cells along R (number of nodes=nr+1)"
    )
    nz: int | None = Field(
        default=None, description="Number of cells along Z (number of nodes=nz+1)"
    )
    n_azimuthal_modes: int | None = Field(
        default=None, description="Number of azimuthal modes"
    )
    rmin: float | None = Field(
        default=None, description="Position of first node along R [m]"
    )
    rmax: float | None = Field(
        default=None, description="Position of last node along R [m]"
    )
    zmin: float | None = Field(
        default=None, description="Position of first node along Z [m]"
    )
    zmax: float | None = Field(
        default=None, description="Position of last node along Z [m]"
    )
    bc_rmin: str | None = Field(
        default=None,
        description="Boundary condition at min R: One of open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_rmax: str | None = Field(
        default=None,
        description="Boundary condition at max R: One of open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_zmin: str | None = Field(
        default=None,
        description="Boundary condition at min Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_zmax: str | None = Field(
        default=None,
        description="Boundary condition at max Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    moving_window_velocity: list[float] | None = Field(
        default=None, description="Moving frame velocity [m/s]"
    )
    refined_regions: list = Field(
        default_factory=list,
        description="List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor], with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc, lo and hi being vectors of length 2 specifying the extent of the region, and refinement_factor defaulting to [2,2] (relative to next lower level)",
    )
    lower_bound_particles: list[float] | None = Field(
        default=None, description="Position of particle lower bound [m]"
    )
    upper_bound_particles: list[float] | None = Field(
        default=None, description="Position of particle upper bound [m]"
    )
    rmin_particles: float | None = Field(
        default=None, description="Position of min particle boundary along R [m]"
    )
    rmax_particles: float | None = Field(
        default=None, description="Position of max particle boundary along R [m]"
    )
    zmin_particles: float | None = Field(
        default=None, description="Position of min particle boundary along Z [m]"
    )
    zmax_particles: float | None = Field(
        default=None, description="Position of max particle boundary along Z [m]"
    )
    # --Like bc_rmin, the radial entry may be None since the lower radial boundary will usually be the axis.
    lower_boundary_conditions_particles: list[str | None] | None = Field(
        default=None,
        description="Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal",
    )
    upper_boundary_conditions_particles: list[str | None] | None = Field(
        default=None,
        description="Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal",
    )
    bc_rmin_particles: str | None = Field(
        default=None,
        description="Boundary condition at min R for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_rmax_particles: str | None = Field(
        default=None,
        description="Boundary condition at max R for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_zmin_particles: str | None = Field(
        default=None,
        description="Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_zmax_particles: str | None = Field(
        default=None,
        description="Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal",
    )
    guard_cells: list[int] | None = Field(
        default=None, description="Number of guard cells used along each direction"
    )
    pml_cells: list[int] | None = Field(
        default=None,
        description="Number of Perfectly Matched Layer (PML) cells along each direction",
    )

    @model_validator(mode="after")
    @resolve_once
    def _resolve_grid(self) -> Self:
        # Sanity check and init of input arguments related to grid parameters
        assert (self.number_of_cells is not None) or (
            self.nr is not None and self.nz is not None
        ), "Either number_of_cells or nr and nz must be specified"
        assert (self.lower_bound is not None) or (
            self.rmin is not None and self.zmin is not None
        ), "Either lower_bound or rmin and zmin must be specified"
        assert (self.upper_bound is not None) or (
            self.rmax is not None and self.zmax is not None
        ), "Either upper_bound or rmax and zmax must be specified"
        # --Allow bc_rmin to be None since it will usually be the axis.
        assert (self.lower_boundary_conditions is not None) or (
            self.bc_zmin is not None
        ), "Either lower_boundary_conditions or bc_rmin and bc_zmin must be specified"
        assert (self.upper_boundary_conditions is not None) or (
            self.bc_rmax is not None and self.bc_zmax is not None
        ), "Either upper_boundary_conditions or bc_rmax and bc_zmax must be specified"

        # Resolve and synchronize the vector and per-axis forms, see _PICMIGrid
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        self._resolve_axis_groups()

        # Sanity check on dimensionality of vector quantities
        assert len(self.number_of_cells) == 2, "Wrong number of cells specified"
        assert len(self.lower_bound) == 2, "Wrong number of lower bounds specified"
        assert len(self.upper_bound) == 2, "Wrong number of upper bounds specified"
        assert len(self.lower_boundary_conditions) == 2, (
            "Wrong number of lower boundary conditions specified"
        )
        assert len(self.upper_boundary_conditions) == 2, (
            "Wrong number of upper boundary conditions specified"
        )
        assert len(self.lower_bound_particles) == 2, (
            "Wrong number of particle lower bounds specified"
        )
        assert len(self.upper_bound_particles) == 2, (
            "Wrong number of particle upper bounds specified"
        )
        assert len(self.lower_boundary_conditions_particles) == 2, (
            "Wrong number of lower particle boundary conditions specified"
        )
        assert len(self.upper_boundary_conditions_particles) == 2, (
            "Wrong number of upper particle boundary conditions specified"
        )

        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2, 2])
            assert len(region[1]) == 2, (
                "The lo extent of the refined region must be a vector of length 2"
            )
            assert len(region[2]) == 2, (
                "The hi extent of the refined region must be a vector of length 2"
            )
            assert len(region[3]) == 2, (
                "The refinement factor of the refined region must be a vector of length 2"
            )

        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2, 2]):
        """Add a refined region.

        Parameters
        ----------
        level : int
            The refinement level, with 1 being the first level of refinement, 2 being the second etc.

        lo, hi : list of float
            The lower (lo) and upper (hi) corner of the region, each a vector of length 2

        refinement_factor : list of int, optional
            The refinement factor relative to the next lower level (default: [2, 2])
        """
        # assign (instead of appending in place), so that the region is validated and marked as set
        self.refined_regions = [*self.refined_regions, [level, lo, hi, refinement_factor]]


class PICMI_Cartesian2DGrid(_PICMIGrid):
    """
    Two dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

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
    # vector of values (such as number_of_cells). Both forms are kept in sync (see _PICMIGrid), but the
    # implementation should use the vectors to access the user input.

    number_of_dimensions: ClassVar[int] = 2
    _axis_groups: ClassVar[tuple[_AxisGroup, ...]] = _grid_axis_groups(("x", "y"))

    # Vector forms (the internally-used representation)
    number_of_cells: list[int] | None = Field(
        default=None,
        description="Number of cells along each axis (number of nodes is number_of_cells+1)",
    )
    lower_bound: list[float] | None = Field(
        default=None, description="Position of the node at the lower bound [m]"
    )
    upper_bound: list[float] | None = Field(
        default=None, description="Position of the node at the upper bound [m]"
    )
    lower_boundary_conditions: list[str | None] | None = Field(
        default=None,
        description="Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    upper_boundary_conditions: list[str | None] | None = Field(
        default=None,
        description="Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    # Per-axis scalar forms (resolved into the vector forms during validation)
    nx: int | None = Field(
        default=None, description="Number of cells along X (number of nodes=nx+1)"
    )
    ny: int | None = Field(
        default=None, description="Number of cells along Y (number of nodes=ny+1)"
    )
    xmin: float | None = Field(
        default=None, description="Position of first node along X [m]"
    )
    xmax: float | None = Field(
        default=None, description="Position of last node along X [m]"
    )
    ymin: float | None = Field(
        default=None, description="Position of first node along Y [m]"
    )
    ymax: float | None = Field(
        default=None, description="Position of last node along Y [m]"
    )
    bc_xmin: str | None = Field(
        default=None,
        description="Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_xmax: str | None = Field(
        default=None,
        description="Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_ymin: str | None = Field(
        default=None,
        description="Boundary condition at min Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_ymax: str | None = Field(
        default=None,
        description="Boundary condition at max Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    moving_window_velocity: list[float] | None = Field(
        default=None, description="Moving frame velocity [m/s]"
    )
    refined_regions: list = Field(
        default_factory=list,
        description="List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor], with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc, lo and hi being vectors of length 2 specifying the extent of the region, and refinement_factor defaulting to [2,2] (relative to next lower level)",
    )
    lower_bound_particles: list[float] | None = Field(
        default=None, description="Position of particle lower bound [m]"
    )
    upper_bound_particles: list[float] | None = Field(
        default=None, description="Position of particle upper bound [m]"
    )
    xmin_particles: float | None = Field(
        default=None, description="Position of min particle boundary along X [m]"
    )
    xmax_particles: float | None = Field(
        default=None, description="Position of max particle boundary along X [m]"
    )
    ymin_particles: float | None = Field(
        default=None, description="Position of min particle boundary along Y [m]"
    )
    ymax_particles: float | None = Field(
        default=None, description="Position of max particle boundary along Y [m]"
    )
    lower_boundary_conditions_particles: list[str | None] | None = Field(
        default=None,
        description="Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal",
    )
    upper_boundary_conditions_particles: list[str | None] | None = Field(
        default=None,
        description="Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal",
    )
    bc_xmin_particles: str | None = Field(
        default=None,
        description="Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_xmax_particles: str | None = Field(
        default=None,
        description="Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_ymin_particles: str | None = Field(
        default=None,
        description="Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_ymax_particles: str | None = Field(
        default=None,
        description="Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal",
    )
    guard_cells: list[int] | None = Field(
        default=None, description="Number of guard cells used along each direction"
    )
    pml_cells: list[int] | None = Field(
        default=None,
        description="Number of Perfectly Matched Layer (PML) cells along each direction",
    )

    @model_validator(mode="after")
    @resolve_once
    def _resolve_grid(self) -> Self:
        self._resolve_cartesian_grid()
        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2, 2]):
        """Add a refined region.

        Parameters
        ----------
        level : int
            The refinement level, with 1 being the first level of refinement, 2 being the second etc.

        lo, hi : list of float
            The lower (lo) and upper (hi) corner of the region, each a vector of length 2

        refinement_factor : list of int, optional
            The refinement factor relative to the next lower level (default: [2, 2])
        """
        # assign (instead of appending in place), so that the region is validated and marked as set
        self.refined_regions = [*self.refined_regions, [level, lo, hi, refinement_factor]]


class PICMI_Cartesian3DGrid(_PICMIGrid):
    """
    Three dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

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
    # vector of values (such as number_of_cells). Both forms are kept in sync (see _PICMIGrid), but the
    # implementation should use the vectors to access the user input.

    number_of_dimensions: ClassVar[int] = 3
    _axis_groups: ClassVar[tuple[_AxisGroup, ...]] = _grid_axis_groups(("x", "y", "z"))

    # Vector forms (the internally-used representation)
    number_of_cells: list[int] | None = Field(
        default=None,
        description="Number of cells along each axis (number of nodes is number_of_cells+1)",
    )
    lower_bound: list[float] | None = Field(
        default=None, description="Position of the node at the lower bound [m]"
    )
    upper_bound: list[float] | None = Field(
        default=None, description="Position of the node at the upper bound [m]"
    )
    lower_boundary_conditions: list[str | None] | None = Field(
        default=None,
        description="Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    upper_boundary_conditions: list[str | None] | None = Field(
        default=None,
        description="Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    # Per-axis scalar forms (resolved into the vector forms during validation)
    nx: int | None = Field(
        default=None, description="Number of cells along X (number of nodes=nx+1)"
    )
    ny: int | None = Field(
        default=None, description="Number of cells along Y (number of nodes=ny+1)"
    )
    nz: int | None = Field(
        default=None, description="Number of cells along Z (number of nodes=nz+1)"
    )
    xmin: float | None = Field(
        default=None, description="Position of first node along X [m]"
    )
    xmax: float | None = Field(
        default=None, description="Position of last node along X [m]"
    )
    ymin: float | None = Field(
        default=None, description="Position of first node along Y [m]"
    )
    ymax: float | None = Field(
        default=None, description="Position of last node along Y [m]"
    )
    zmin: float | None = Field(
        default=None, description="Position of first node along Z [m]"
    )
    zmax: float | None = Field(
        default=None, description="Position of last node along Z [m]"
    )
    bc_xmin: str | None = Field(
        default=None,
        description="Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_xmax: str | None = Field(
        default=None,
        description="Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_ymin: str | None = Field(
        default=None,
        description="Boundary condition at min Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_ymax: str | None = Field(
        default=None,
        description="Boundary condition at max Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_zmin: str | None = Field(
        default=None,
        description="Boundary condition at min Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    bc_zmax: str | None = Field(
        default=None,
        description="Boundary condition at max Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann",
    )
    moving_window_velocity: list[float] | None = Field(
        default=None, description="Moving frame velocity [m/s]"
    )
    refined_regions: list = Field(
        default_factory=list,
        description="List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor], with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc, lo and hi being vectors of length 3 specifying the extent of the region, and refinement_factor defaulting to [2,2,2] (relative to next lower level)",
    )
    lower_bound_particles: list[float] | None = Field(
        default=None, description="Position of particle lower bound [m]"
    )
    upper_bound_particles: list[float] | None = Field(
        default=None, description="Position of particle upper bound [m]"
    )
    xmin_particles: float | None = Field(
        default=None, description="Position of min particle boundary along X [m]"
    )
    xmax_particles: float | None = Field(
        default=None, description="Position of max particle boundary along X [m]"
    )
    ymin_particles: float | None = Field(
        default=None, description="Position of min particle boundary along Y [m]"
    )
    ymax_particles: float | None = Field(
        default=None, description="Position of max particle boundary along Y [m]"
    )
    zmin_particles: float | None = Field(
        default=None, description="Position of min particle boundary along Z [m]"
    )
    zmax_particles: float | None = Field(
        default=None, description="Position of max particle boundary along Z [m]"
    )
    lower_boundary_conditions_particles: list[str | None] | None = Field(
        default=None,
        description="Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal",
    )
    upper_boundary_conditions_particles: list[str | None] | None = Field(
        default=None,
        description="Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal",
    )
    bc_xmin_particles: str | None = Field(
        default=None,
        description="Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_xmax_particles: str | None = Field(
        default=None,
        description="Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_ymin_particles: str | None = Field(
        default=None,
        description="Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_ymax_particles: str | None = Field(
        default=None,
        description="Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_zmin_particles: str | None = Field(
        default=None,
        description="Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal",
    )
    bc_zmax_particles: str | None = Field(
        default=None,
        description="Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal",
    )
    guard_cells: list[int] | None = Field(
        default=None, description="Number of guard cells used along each direction"
    )
    pml_cells: list[int] | None = Field(
        default=None,
        description="Number of Perfectly Matched Layer (PML) cells along each direction",
    )

    @model_validator(mode="after")
    @resolve_once
    def _resolve_grid(self) -> Self:
        self._resolve_cartesian_grid()
        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2, 2, 2]):
        """Add a refined region.

        Parameters
        ----------
        level : int
            The refinement level, with 1 being the first level of refinement, 2 being the second etc.

        lo, hi : list of float
            The lower (lo) and upper (hi) corner of the region, each a vector of length 3

        refinement_factor : list of int, optional
            The refinement factor relative to the next lower level (default: [2, 2, 2])
        """
        # assign (instead of appending in place), so that the region is validated and marked as set
        self.refined_regions = [*self.refined_regions, [level, lo, hi, refinement_factor]]


PICMI_AnyGrid = (
    PICMI_CylindricalGrid
    | PICMI_Cartesian1DGrid
    | PICMI_Cartesian2DGrid
    | PICMI_Cartesian3DGrid
    | PICMI_Grid
)

_ElectromagneticSolverMethod = Literal[
    "Yee", "CKC", "Lehe", "PSTD", "PSATD", "GPSTD", "DS", "ECT"
]


class PICMI_ElectromagneticSolver(PICMI_Solver):
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

    # Retained for backwards compatibility reasons.
    # The type annotation of `method` is the ground-truth.
    methods_list: ClassVar[list[str]] = list(get_args(_ElectromagneticSolverMethod))

    grid: PICMI_AnyGrid = Field(description="Grid object for the diagnostic")
    method: _ElectromagneticSolverMethod | None = Field(
        default=None,
        description="The advance method use to solve Maxwell's equations. The default method is code dependent.",
    )
    stencil_order: list[int] | None = Field(
        default=None, description="Order of stencil for each axis (-1=infinite)"
    )
    cfl: float | None = Field(
        default=None, description="Fraction of the Courant-Friedrich-Lewy criteria [1]"
    )
    source_smoother: PICMI_BinomialSmoother | None = Field(
        default=None, description="Smoother object to apply to the sources"
    )
    field_smoother: PICMI_BinomialSmoother | None = Field(
        default=None, description="Smoother object to apply to the fields"
    )
    subcycling: int | None = Field(
        default=None, description="Level of subcycling for the GPSTD solver"
    )
    galilean_velocity: list[float] | None = Field(
        default=None, description="Velocity of Galilean reference frame [m/s]"
    )
    divE_cleaning: bool | None = Field(
        default=None, description="Solver uses div(E) cleaning if True"
    )
    divB_cleaning: bool | None = Field(
        default=None, description="Solver uses div(B) cleaning if True"
    )
    pml_divE_cleaning: bool | None = Field(
        default=None, description="Solver uses div(E) cleaning in the PML if True"
    )
    pml_divB_cleaning: bool | None = Field(
        default=None, description="Solver uses div(B) cleaning in the PML if True"
    )


_ElectrostaticSolverMethod = Literal["FFT", "Multigrid"]


class PICMI_ElectrostaticSolver(PICMI_Solver):
    """
    Electrostatic field solver
    """

    # Retained for backwards compatibility reasons.
    # The type annotation of `method` is the ground-truth.
    methods_list: ClassVar[list[str]] = list(get_args(_ElectrostaticSolverMethod))

    grid: PICMI_AnyGrid = Field(description="Grid object for the diagnostic")
    method: _ElectrostaticSolverMethod | None = Field(
        default=None, description="One of 'FFT', or 'Multigrid'"
    )
    required_precision: float | None = Field(
        default=None, description="Level of precision required for iterative solvers"
    )
    maximum_iterations: int | None = Field(
        default=None, description="Maximum number of iterations for iterative solvers"
    )


_MagnetostaticSolverMethod = Literal["FFT", "Multigrid"]


class PICMI_MagnetostaticSolver(PICMI_Solver):
    """
    Magnetostatic field solver
    """

    # Retained for backwards compatibility reasons.
    # The type annotation of `method` is the ground-truth.
    methods_list: ClassVar[list[str]] = list(get_args(_MagnetostaticSolverMethod))

    grid: PICMI_AnyGrid = Field(description="Grid object for the diagnostic")
    method: _MagnetostaticSolverMethod | None = Field(
        default=None, description="One of 'FFT', or 'Multigrid'"
    )


PICMI_AnySolver = (
    PICMI_ElectromagneticSolver
    | PICMI_ElectrostaticSolver
    | PICMI_MagnetostaticSolver
    | PICMI_Solver
)
