The Accelerators Modeling Interface (AMI) Standard
==================================================

VERSION: **1.0.0** (October 12th, 2016)

Introduction
------------

The AMI standard establishes conventions for the naming and structuring of input files for particle accelerator simulations. 

There is already broad acceptance of the Standard Input Format, also known as the MAD format that was developed for the MAD and MAD-X codes. Given this wide acceptance, the AMI standard incorporates many of the MAD conventions to describe beamline elements but also expands upon them significantly to include new beamline elements and new functionalities. For example, the MAD format utilizes a quadrupole coefficient for an idealized quadrupole. AMI subsumes this but also includes new functionality to treat quadrupoles with realistic fringe and include nonlinear effects associated with high-order terms in the vector potential. This is accomplished by defining a new attribute, gradients, that points to a file that contains tabulated values of the generalized gradients. This allows the computation of nonlinear maps with realistic field profiles. Also, it allows for inclusion of nonlinear effects (other than kinematic nonlinearities) when doing particle tracking.

Note that syntactic issues such as whether or not an input record contains a semicolon (which is one of the features that distinguishes MAD-X from MAD) are not part of AMI, as that is considered an implementation issue.

A list of beamline elements follows. Note that, while AMI includes many extensions to the MAD format, there is certain new functionality that is common to many elements. These extensions are:

- some thick elements (dipoles, quadrupoles,...) have an optional "gradients" parameter that enables input of the generalized gradients associated with realistic models of the elements.
- all thick elements have an optional "slices" parameter. This is typically used for diagnostic purposes, e.g., to print a diagnostic quantity (rms size, rms emittance, etc.) at the beginning or middle or end of every slice. But it can also be used for other purposes. If this parameter is omitted then the default number of slices is 1.
- all thick elements have an optional "method" parameter that describes how the element is to be used by the code. For example, "method" can describe the algorithm that will be used to track particles through the element. Or, "method" can be used to describe that the code will produce a Lie map or Taylor map for the element. Note that beam dynamics codes normally have a global means for describing this behavior. Hence this parameter is optional, but it provides the user with added flexibility if desired. 

Lattice elements
----------------
  - `Element` - "Fundamental element of accelerator lattice."
    - **type**: *object*
    - `Length=length=L=l` - **type**: *double* - "Length of element [m]."
    - `Offsetx=offsetx=Ox=ox` - **type**: *double* - "Offset along X [m]."
    - `Offsety=offsety=Oy=oy` - **type**: *double* - "Offset along Y [m]."
    - `Offsetz=offsetz=Oz=oz` - **type**: *double* - "Offset along Z [m]."
    - `Tiltx=tiltx=Tx=tx` - **type**: *double* - "Tilt along X [degree]."
    - `Tilty=tilty=Tx=ty` - **type**: *double* - "Tilt along Y [degree]."
    - `Tiltz=tiltz=Tx=tz` - **type**: *double* - "Tilt along Z [degree]."
    - `Aperture=aperture=Rmax=rmax` - **type**: *double* - "Aperture [m]."
    
  - `Marker` - "Marker along the lattice."  
    - **type**: *object*
    
  - `Drift` - "Drift."
    - **type**: *Element*
    
  - `Solen` - "Magnetic Solenoid."
    - **type**: *Element*
    - `Bz` - **type**: *double* - "Bz field [V/m]."

  - `Dipo` - "Dipole (magnetic or electric)."
    - **type**: *Element*
    - `Ex` - **type**: *double* - "Ex field [V/m]."
    - `Ey` - **type**: *double* - "Ey field [V/m]."
    - `Bx` - **type**: *double* - "Bx field [T]."
    - `By` - **type**: *double* - "By field [T]."
    
  - `Quad` - "Quadrupole (magnetic or electric)."
    - **type**: *Element*
    - `DE` - **type**: *double* - "Electric field gradient [V/m**2]."
    - `DB` - **type**: *double* - "Magnetic field gradient [T/m]."
    
  - `Sext` - "Sextupole (magnetic or electric)."
    - **type**: *Element*
    - `DE` - **type**: *double* - "Electric field gradient [V/m**2]."
    - `DB` - **type**: *double* - "Magnetic field gradient [T/m]."

  - `RFC` - "Radio Frequency Cavity"
    - **type**: *Element*
    - `freq` - **type**: *double* - "Frequency [Hz]."
    - `phase` - **type**: *double* - "Nominal phase [degree]."
    - `Vmax` - **type**: *double* - "Maximum integrated voltage [V]."
    
  - `Undulator` - "Magnetic Undulator"
    - **type**: *Element*
    - `period` - **type**: *double* - "Period [m]."
    - `Bxmax` - **type**: *double* - "Maximum on-axis Bx field [T]."
    - `Bxphase` - **type**: *double* - "Horizontal phase [degree]."
    - `Bymax` - **type**: *double* - "Maximum on-axis By field [T]."
    - `Byphase` - **type**: *double* - "Vertical phase [degree]."
