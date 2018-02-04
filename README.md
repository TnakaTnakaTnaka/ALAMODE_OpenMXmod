
ALAMODE_OpenMX_tools
====

# Overview

ALAMODE is a software package to analyze phonon properties by using external DFT calculation code, OpenMX (DFT calcuation code) is not supported to ALAMODE. disp_openmx_mod.py and extract_openmx_mod.py are script to use ALAMODE by using OpenMX code.

supercell.py is a script to generate supercell data file of OpenMX.


# Tools
## Original scripts and module script

* displace.py - script to create DFT calculation input file for all displacement pattern in ALAMODE. This script is developed by ALAMODE developer, a part of code is add by Yuto Tanaka to support OpenMX.

* extract.py - script to extract displacement, force, energy from output files of DFT calculations. This script is developed by ALAMODE developer, a part of code is add by Yuto Tanaka to support OpenMX.

* disp_openmx_mod.py - additional module of displace.py to support OpenMX in ALAMODE.

* extract_openmx_mod.py - additional module fo extract.py to support OpenMX in ALAMODE.

## Supercell script
* supercell.py - Simple script to generate supercell data file. OpenMX is supported.


# Reference
ALAMODE software package

http://alamode.readthedocs.io/


# Author
Yuto Tanaka (Kanazawa University)


