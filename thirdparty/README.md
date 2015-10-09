LEMOS Rostock Extensions
====================
Collections of OpenFOAM-Extensions created at the Institute for
Modeling and Numerical Simulation at the University of Rostock

The thirparty directory contains all LEMOS developments which 
depend on external or thirdparty software.
If the corresponding dependencies are added to the original OpenFOAM 
framework someday, the depended tools, libraries or applications, created 
by LEMOS, will be transfered to the standard LEMOS libraries and 
application folders.

Due to external dependencies the thirdparty is not build automatically.
Each application, tool or library must be build manually, executing the 
Allwmake script in the corresponding folder. The tools are not part of
the libLEMOS library, meaning that their own libraries have to be included.

Installation
============

1. If not already done source  ". $FOAM_SRC/LEMOS-2.4.x/bashrc"
2. Execute the Allwmake in the correponding application folder to build the 
   libraries and/or executables.
