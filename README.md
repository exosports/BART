Transit
=======

A radiative-transfer code for planetary atmospheres

This code was orginially written by Patricio Rojo in his PhD thesis with Joseph Harrington at Cornell University.  It only handled the case of tangent-geometry absorption at the limb, as is the case for exoplanet transit observations, but it is generally written in an object-oriented style in C, as we planned to extend it as needed for other situations.  Harrington's current group at the University of Central Florida is now modifying it to handle emission (as in exoplanet secondary eclipses) as well.  It will be used as a component in the group's Bayesian Atmospheric Radiative Transfer (BART) project.

The new version of transit is being developed by:

Joseph Harrington  
Jasmina Blecic  
Patricio Cubillos  
Andrew Foster  

With support from:

Jonathan Fortney (UCSC), Nikku Madhusudhan (Yale), Patricio Rojo (U. de Chile)

Test after each upgrade, record output, compare output of next upgrade
with upgrade turned off to previous output to make sure we didn't
wreck anything.
