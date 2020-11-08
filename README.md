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
* [Patricio Cubillos](https://github.com/pcubillos/) (UCF, IWF) <patricio.cubillos@oeaw.ac.at>
* Jasmina Blecic (UCF)
* Joseph Harrington (UCF)
* Patricio Rojo (U. de Chile)
* Nate Lust (UCF)
* Oliver Bowman (UCF)
* Madison Stemm (UCF)
* Andrew Foster (UCF)
* Ryan Challener (UCF)
* Michael Himes (UCF)

With support from:
* Thomas J. Loredo (Cornell)
* Jonathan Fortney (UCSC)
* Nikku Madhusudhan (Yale)

### Getting Started:
You can find the BART User Manual [here](https://exosports.github.io/BART/doc/BART_User_Manual.html). If you run into trouble while using BART, please feel free to use our [BART User Mailing List](https://physics.ucf.edu/mailman/listinfo/bart-user). You'll have to sign up to post (for spam reasons), but we encourage you to join. If you're a developer interested in contributing to BART, you may also be interested in our [BART Developer Mailing List](https://physics.ucf.edu/mailman/listinfo/bart-devel) (see CONTRIBUTING for more information).

### Install and Compile:
Download the latest stable version from the BART
[releases](https://github.com/exosports/BART/releases) page
(TBD).  Alternatively, clone the repository to your local
machine with the following terminal commands.
First create a working directory to place the code:
```shell
mkdir BART_demo/
cd BART_demo/
topdir=`pwd`
```

Clone the repository with all its submodules:
```shell
git clone --recursive https://github.com/exosports/BART $topdir/BART/
```

Enter the BART directory and build the conda environment:
```shell
cd $topdir/BART/
conda env create -f environment.yml
conda activate bart
```

Compile the transit module programs:
```shell
cd $topdir/BART/modules/transit/
make
```

Compile the MCcubed routines:
```shell
cd $topdir/BART/modules/MCcubed/
make
```

To remove the program binaries, execute (in the respective
directories):
```shell
make clean
```

### Quick Example:

The following script lets you quickly fit a methane emission spectrum model to a set of 10 filters between 2 and 4 um.  These instructions are meant to be executed from the shell terminal.  To begin, follow the instructions in the previous Section to install and compile the code.  Now, create a working directory in your top directory to place the files and execute the programs:
```shell
cd $topdir
mkdir run/
cd run/
```

Download the methane line-transition database from the HITRAN server:
```shell
wget --user=HITRAN --password=getdata -N https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/06_hit12.zip 
unzip 06_hit12.zip
```

Copy the pylineread configuration file and run pylineread to make the transition-line-information (TLI) file:
```shell
cp $topdir/BART/examples/demo/pyline_demo.cfg .  
$topdir/BART/modules/transit/pylineread/src/pylineread.py -c pyline_demo.cfg
```

Copy the transit configuration file and run it to make a table of opacities:
```shell
cp $topdir/BART/examples/demo/transit_demo.cfg .  
$topdir/BART/modules/transit/transit/transit -c transit_demo.cfg --justOpacity
```

Copy and run the BART configuration file for eclipse geometry:
```shell
cp $topdir/BART/examples/demo/BART_eclipse.cfg .
$topdir/BART/BART.py -c BART_eclipse.cfg
```

Copy and run the BART configuration file for transit geometry:
```shell
cp $topdir/BART/examples/demo/BART_transit.cfg .
$topdir/BART/BART.py -c BART_transit.cfg
```


### Be Kind:

A number of graduate and undergraduate-student dreams and hopes have been sacrificed on the making of this project.  Please, be kind and aknowledge their effort by citing the articles asociated to this project:

  [Cubillos et al. 2016: The Bayesian Atmospheric Radiative-Transifer Code for Exoplanet Modeling](), in preparation.  
  [Blecic et al. 2016: The Bayesian Atmospheric Radiative-Transifer Code for Exoplanet Modeling](), in preparation.  
  [Harrington et al. 2016: The Bayesian Atmospheric Radiative-Transifer Code for Exoplanet Modeling](), in preparation.  

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

Copyright (C) 2015-2016 University of Central Florida.
All rights reserved.  

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
Patricio Cubillos <patricio.cubillos@oeaw.ac.at>  
Jasmina Blecic <jasmina@physics.ucf.edu>  
Joseph Harrington <jh@physics.ucf.edu>  

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

<!--
**STATUS OF COMPONENTS**

| Component     | Code          | Doc           | Package | Git | Centralgit(here) |
| ------------- | --------------| ------------  |---------|-----| -----------------|
| comm loop     | done          |               |         | yes |                  |
| MC-cubed      | done          | To be revised |         | yes |                  |
| TEA           | done          | On revision   |         | yes |                  |
| transit       | done          |               |         | yes |                  |
| BART          | done          |               |         | yes |                  |
-->
