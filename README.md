# gromacs_AFM_FH-MD

The code simulates the Atomic Force Microscopy (AFM) via Fluctuating Hydrodynamics-Molecular dynamics (FH/MD) method. To run this code, the input file, initial configuration file and velocity file are needed. In the input file, the number of FH grids and alpha and beta should be specified. The initial configuration file should contain a substrate, a spherical tip and water solution and corresponding forcefield. The velocity of tip is specified in the velocity file. The output file contains instantaneous force on the tip.

## Installation

The code can be installed as original GROMACS 2016.4, see [INSTALL](https://github.com/ikorotkin/gromacs_AFM_FH-MD/blob/main/INSTALL) file for details.

## Licensing

-   GROMACS is open source software distributed under the GNU Lesser General Public License (LGPL) Version 2.1 or (at your option) any later version. GROMACS includes optional code covered by several different licences. See [COPYING](https://github.com/ikorotkin/gromacs_AFM_FH-MD/blob/main/COPYING) file for details.
-   AFM FH/MD add-on to GROMACS is full open source under [MIT license](https://github.com/ikorotkin/gromacs_AFM_FH-MD/blob/main/LICENSE_FH-MD).
