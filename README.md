### BART
> Bayesian Atmospheric Radiative Transfer fitting code

This code implements a Bayesian, Monte Carlo-driven,
radiative-transfer scheme for extracting parameters from spectra of
planetary atmospheres.  

### Table of Contents:
* [Team Members](#team-members)
* [Getting Started](#getting-started)
* [Install and Compile](#install-and-compile)
* [Quick Example](#quick-example)
* [Be Kind](#be-kind)
* [License](#license)


### Team Members:
* [Patricio Cubillos](https://github.com/pcubillos/) (UCF) <pcubillos@fulbrightmail.org>
* Jasmina Blecic (UCF)
* Joseph Harrington (UCF)
* Patricio Rojo (U. de Chile)
* Oliver Bowman (UCF)
* Madison Stemm (UCF)
* Andrew Foster (UCF)

With support from:
* Thomas J. Loredo (Cornell)
* Jonathan Fortney (UCSC)
* Nikku Madhusudhan (Yale)

### Getting Started:
Get the BART user's manual [here](doc/) (TBD).

### Install and Compile:
Download the latest stable version from the BART
[releases](https://github.com/joeharr4/BART/releases) page
(TBD).  Alternatively, clone the repository to your local
machine with the following terminal commands.
First create a working directory to place the code:
```shell
cd
mkdir tmp/
mkdir tmp/BART_demo/
cd tmp/BART_demo/
```

Clone the repository with all its submodules:
```shell
git clone --recursive https://github.com/joeharr4/BART BART/
```

Compile the transit module programs:
```shell
barttop=`pwd`
cd BART/modules/transit/pylineread/src/fortran/
make
cd ../../../pu/
make
cd ../transit/
make
./config
./compile
```

Compile the MCcubed routines:
```shell
cd ../../MCcubed/src/cfuncs/
make
```

To remove the program binaries, execute (in the respective
directories):
```shell
make clean
```

### Quick Example:

### Be Kind:

### License:

Bayesian Atmospheric Radiative Transfer (BART), a code to infer
properties of planetary atmospheres based on observed spectroscopic
information.  

This project was completed with the support of the NASA Planetary
Atmospheres Program, grant NNX12AI69G, held by Principal Investigator
Joseph Harrington. Principal developers included graduate students
Patricio E. Cubillos and Jasmina Blecic, programmer Madison Stemm, and
undergraduates M. Oliver Bowman and Andrew S. D. Foster.  The included
'transit' radiative transfer code is based on an earlier program of
the same name written by Patricio Rojo (Univ. de Chile, Santiago) when
he was a graduate student at Cornell University under Joseph
Harrington.  Statistical advice came from Thomas J. Loredo and Nate
B. Lust.  

Copyright (C) 2015 University of Central Florida.  All rights reserved.  

This is a test version only, and may not be redistributed to any third
party.  Please refer such requests to us.  This program is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  

Our intent is to release this software under an open-source,
reproducible-research license, once the code is mature and the first
research paper describing the code has been accepted for publication
in a peer-reviewed journal.  We are committed to development in the
open, and have posted this code on github.com so that others can test
it and give us feedback.  However, until its first publication and
first stable release, we do not permit others to redistribute the code
in either original or modified form, nor to publish work based in
whole or in part on the output of this code.  By downloading, running,
or modifying this code, you agree to these conditions.  We do
encourage sharing any modifications with us and discussing them
openly.  

We welcome your feedback, but do not guarantee support.  Please send
feedback or inquiries to:  
Patricio Cubillos <pcubillos[at]fulbrightmail.org>  
Jasmina Blecic <jasmina[at]physics.ucf.edu>  
Joseph Harrington <jh[at]physics.ucf.edu>  

or alternatively,  
Joseph Harrington, Patricio Cubillos, and Jasmina Blecic  
UCF PSB 441  
4111 Libra Drive  
Orlando, FL 32816-2385  
USA  

Thank you for testing BART!  

### BART Development Plan

0. Fully understand Pato's code (all, esp. Patricio, Jasmina) (DONE!)
1. Split transit code into initialization and loop (Patricio, Andrew) (DONE!)
2. Drive loop with arbitrary parameters (Patricio)
3. Interface transit code with MC driver (Jasmina, Patricio) (READY for testing)
4. First test: H2O-only, transit, constant abundances, constant T(p)
   (good enough for NESSF!  15 Jan)
5. Upgrades required to publish (do in any order):
   many molecules' line lists (Patricio)
   eclipse geometry (Jasmina) (READY for tests)
   CEA initial conditions/constant scaling (Jasmina) (READY for tests)
   Madhu T(p) (Jasmina)
6. Validation vs. Madhu
   (can do science now)
7. Performance upgrades:
   correlated-k
   multiple chains from one node driving single RT code
   multiple chains from one node driving multiple RT codes
   multiple chains from multiple nodes driving multiple RT codes
8. Science upgrades:
   flux-balanced layers (a la Burrows)
   constant abundances/T(p) throughout atmosphere
   line-list-difference tests
   dayside map
   phase curve
   dayside map + phase curve

Test after each upgrade, record output, compare output of next upgrade
with upgrade turned off to previous output to make sure we didn't
wreck anything.

**STATUS OF COMPONENTS**

| Component     | Code          | Doc           | Package | Git | Centralgit(here) |
| ------------- | --------------| ------------  |---------|-----| -----------------|
| comm loop     | done          |               |         | yes |                  |
| MC-cubed      | done          | To be revised |         | yes |                  |
| TEA           | done          | On revision   |         | yes |                  |
| transit       | dev           |               |         | yes |                  |
| BART          | dev           |               |         | dev |                  |

