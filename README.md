# Quick Start
 
# Disclaimer and Licensing
 
yasw is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
yasw is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with yasw. If not, see <http://www.gnu.org/licenses/>.
The full software license for yasw version 1.0.0 
can be found in
file COPYING. 

# Please cite this work

yasw is a free and open source implementation of the 
DCA algorithm for models of strongly correlated electrons. 
The full software license for yasw version 1.0.0 
can be found in
file COPYING. 
You are welcomed to use it and publish data 
obtained with yasw. If you do, please cite this
work. Explain How To Cite This Work. FIXME. TBW.


# Hash of the latest commit 

Hash of the latest commit is also posted at
https://web.ornl.gov/~gz1/hashes.html

# Building and Running yasw

## Required Software

* GNU C++
* The LAPACK and BLAS libraries
* The GSL library
* PsimagLite (see below)

## Optional Software

* make or gmake (only needed to use the Makefile)
* perl (may be needed to run some auxiliary script) 

## Quick Start

1. Use your distribution repository tool to install gcc with support for C++,
the LAPACK and BLAS libraries, the gsl library, make, perl, doxygen and git 
if you don't have them.

2. Issue

    cd someDirectory/

    git clone https://github.com/g1257/PsimagLite.git

    git clone https://github.com/g1257/yasw.git

3. Compile PsimagLite

    cd PsimagLite/lib/

    make -f Makefile.sample

    cd ../../

4. Now issue

    cd yasw

    make

5. You can run it with

   FIXME


