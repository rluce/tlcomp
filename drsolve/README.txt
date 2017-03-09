  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT.                                                        *
  ***************************************************************************

   AUTHORS:

       Antonio Arico', Giuseppe Rodriguez
       University of Cagliari, Italy
       E-mail: arico@unica.it, rodriguez@unica.it

   REFERENCE:

    -  A fast solver for linear systems with displacement structure
       NUMERICAL ALGORITHMS, 55 (2010), pp. 529-556

   SOFTWARE REVISION DATE:

       V1.0, April 2010

   SOFTWARE LANGUAGE:

       MATLAB 7.9 (R2009b) and C

======================================================================
SOFTWARE
======================================================================

This package provides a set of Matlab functions to solve linear
systems whose matrix belongs to the following classes: Cauchy-like,
Toeplitz-like, Toeplitz+Hankel-like and Vandermonde-like. The solvers
have a computational cost O(rn^2), and use O(rn) memory locations,
where n is the size and r is the displacement rank of the matrix. The
Cauchy-like solver applies the generalized Schur algorithm to a
suitable augmented Cauchy-like matrix. The other solvers perform fast
transformations to convert the input system to a Cauchy-like one, and
then apply the Cauchy-like solver. Various pivoting techniques are
implemented.
See the help page of each function for the calling syntax.


======================================================================
PACKAGE
======================================================================

The main directory contains the following files

README.txt     :  This file
Contents.m     :  Contents file
clsolve.m      :  Solution of a Cauchy-like linear system
tsolve.m       :  Solution of a Toeplitz linear system
tlsolve.m      :  Solution of a Toeplitz-like linear system
thsolve.m      :  Solution of a Toeplitz+Hankel linear system
thlsolve.m     :  Solution of a Toeplitz+Hankel-like linear system
vsolve.m       :  Solution of a Vandermonde linear system
vlsolve.m      :  Solution of a Vandermonde-like linear system
t2cl.m         :  Converts a Toeplitz matrix to Cauchy-like
t2tl.m         :  Converts a Toeplitz matrix to Toeplitz-like
tl2cl.m        :  Converts a Toeplitz-like matrix to Cauchy-like
th2cl.m        :  Converts a square Toeplitz+Hankel matrix to Cauchy-like
thl2cl.m       :  Converts a square Toeplitz+Hankel-like matrix to Cauchy-like
v2cl.m         :  Converts a square Vandermonde matrix to Cauchy-like
vl2cl.m        :  Converts a square Vandermonde-like matrix to Cauchy-like
ftimes.m       :  Fourier matrix product
ctimes.m       :  Cosine matrix product
stimes.m       :  Sine matrix product
ttimes.m       :  Toeplitz matrix product
cltimes.m      :  Cauchy-like matrix product
cl2full.m      :  Construct a Cauchy-like matrix from its generators
nroots1.m      :  n-th roots of unity

test/          :  Directory containing the scripts for numerical experiments
src/           :  Directory containing the C source code and some executables
validate/      :  Directory containing the validation scripts


======================================================================
HOW TO INSTALL
======================================================================

When the archive file containing the package is uncompressed, it
creates a directory named drsolve. This directory must be added to the
Matlab search path, either by the "addpath" command, or using the
menus available in the graphical user interface.

Then, the user must change the current directory to drsolve/validate,
and execute the script "installmex", which will install the
appropriate MEX version of clsolve, if available, in the main drsolve
directory. If the MEX file is not available, the Matlab version of
clsolve will be used for computation. The MEX file can be compiled by
the user; see the "HOW TO COMPILE" section.

To test the installation, run the "validate" script in the
drsolve/validate directory. It will check that the installation
process is complete and that the package is working.


======================================================================
HOW TO COMPILE
======================================================================

The source C code of clsolve is contained in the drsolve/src
directory. To compile it one must must copy the files fort.c and
fort.h from the ${MATLABROOT}/extern/examples/refbook directory to
drsolve/src. The string ${MATLABROOT} denotes the directory where
Matlab is installed.

Then, one must issue, in the same directory, the command
	mex clsolve.c fort.c -lmwblas -lmwlapack
and copy the executable in the main drsolve directory.
The MEX program uses the versions of the BLAS and LAPACK libraries
distributed with Matlab.

On UNIX SYSTEMS the process is straightforward, if Matlab and a
suitable C compilers are correctly installed.

We provide a Makefile, if the Unix "make" command is available. In
this case, it is enough to issue one of the commands
	make linux32		(on Linux 32 bit)
	make linux64		(on Linux 64 bit)
	make macosx		(on Mac OS X)
to compile the executable and install it. The Makefile can be easily
tailored to any Unix system.

On WINDOWS SYSTEMS, Matlab officially supports only some commercial
compilers. It is not possible to use the minimal C compiler "lcc"
distributed with Matlab for Windows, as it does not support complex
variables.

Anyway, the MEX version of clsolve can be compiled using the MinGW
porting of the GNU-C compiler (http://www.mingw.org/) and the Gnumex
"MEX configurator" (http://gnumex.sourceforge.net/). In this case it
is necessary to apply a small patch to the gnumex.m file, by modifying
the lines 886-887 from
	mexlibs = {'libmx', 'libmex', 'libmat'};
	englibs = {'libmx', 'libeng', 'libmat'};
to
	mexlibs = {'libmx', 'libmex', 'libmat', 'libmwblas', 'libmwlapack'};
	englibs = {'libmx', 'libeng', 'libmat', 'libmwblas', 'libmwlapack'};
to link the needed libraries.

An IMPORTANT TECHNICAL REMARK which affects the 64-BIT OPERATING
SYSTEMS is the following. Starting from Matlab 7.8, the type of the
integer variables used in the BLAS and LAPACK libraries provided with
Matlab changed from "int" to "ptrdiff_t". This has no effect on 32-bit
architectures, but on 64-bit systems the size of the variables changed
from 4 to 8 bytes. For this reason, to compile clsolve.c, running
Matlab version 7.7 or less on a 64-bit system, it is necessary to
change the line 28 from
	#define PINT	ptrdiff_t	/* starting from Matlab 7.8 (R2009a) */
to
	#define PINT	int		/* up to Matlab 7.7 (R2008b) */


======================================================================
NUMERICAL TESTS
======================================================================

The directory drsolve/test contains some scripts which perform the
numerical experiments described in the paper. 

Some of the tests require the solver described in the paper
   P.C. Hansen and T.F. Chan. Fortran subroutines for general Toeplitz
   systems. ACM Trans. Math. Software, 18(3):256--273, 1992.
A Matlab MEX gateway to this subroutine is available at
   http://bugs.unica.it/~gppe/soft/

