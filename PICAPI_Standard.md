The Particle-In-Cell Application Programming Interface
======================================================

The standard API for controlling particle-in-cell codes

VERSION: **1.0.0** (May 5, 2017)

General access routines
-----------------------------------

  - `get_time()` - **type**: *method*
    - **returns**: time **type**: *double* - "current time of the simulation"

  - `set_time(` - **type**: *method*
    - time - **type** *double* - "updated simulation time"
    - `)`

  - `get_step_number()` - **type**: *method*
    - **returns**: step_number **type**: *integer* - "current step number of the simulation"

  - `set_step_number(` - **type**: *method*
    - step_number - **type** *integer* - "updated simulation step number"
    - `)`

  - `get_step_size()` - **type**: *method*
    - **returns**: step_size **type**: *double* - "current step size of the simulation"

Particle advancing and field interaction routines
-----------------------------------

  - `push_positions(` - **type**: *method*
    - dt - **type** *double* - "time step size"
    - `)`

  - `push_velocities_withE(` - **type**: *method*
    - dt - **type** *double* - "time step size"
    - `)`

  - `push_velocities_withB(` - **type**: *method*
    - dt - **type** *double* - "time step size"
    - `)`

  - `get_self_fields()` - **type**: *method*

  - `get_applied_fields(` - **type**: *method*
    - dt - **type** *double* - "time step size"
    - dtl - **type** *double* - "step size of first half"
    - dtr - **type** *double* - "step size of second half"
    - `)`

  - `calculate_source()` - **type**: *method*

  - `push_Efields(` - **type**: *method*
    - dt - **type** *double* - "time step size"
    - `)`

  - `push_Bfields(` - **type**: *method*
    - dt - **type** *double* - "time step size"
    - `)`

  - `apply_particle_boundary_conditions()` - **type**: *method*

Time step types
------------------

    - `splitting`
      - **type**: *string*
      - **choices**
        - `velocity` - the step will be split so that the position and velocity
                       will be synchronized after a full position advance
         - `position` - the step will be split so that the position and velocity
                        will be synchronized after a full velocity advance

    - `alwayssplit`
      - **type**: *logical*
      - **choices**
        - `True` - "split steps are always done"
        - `False` - "full steps are done for efficiency when possble"

