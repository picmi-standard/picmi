The Accelerators Modeling Interface (AMI) Standard
==================================================

VERSION: **1.0.0** (October 12th, 2016)

Lattice elements
----------------
  - `Element` - "Fundamental element of accelerator lattice."
    - **type**: *object*
    - `Length=length=L=l` - **type**: *double* - "Length of element."
    - `Offsetx=offsetx=Ox=ox` - **type**: *double* - "Offset along X."
    - `Offsety=offsety=Oy=oy` - **type**: *double* - "Offset along Y."
    - `Offsetz=offsetz=Oz=oz` - **type**: *double* - "Offset along Z."
    - `Tiltx=tiltx=Tx=tx` - **type**: *double* - "Tilt along X."
    - `Tilty=tilty=Tx=ty` - **type**: *double* - "Tilt along Y."
    - `Tiltz=tiltz=Tx=tz` - **type**: *double* - "Tilt along Z."
    
  - `Marker` - "Marker along the lattice."  
    - **type**: *object*
    
  - `Drift` - "Drift."
    - **type**: *Element*
    
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
    
