The Accelerators Modeling Interface (AMI) Standard
==================================================

VERSION: **1.0.0** (October 12th, 2016)

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
