BART
====

Bayesian Atmospheric Radiative-Transfer fitting code

This code will implement a Bayesian, Monte Carlo-driven,
radiative-transfer scheme for extracting parameters from spectra of
planetary atmospheres.  It is being developed by:

Joseph Harrington
Jasmina Blecic
Patricio Cubillos
Sarah Blumenthal
Oliver Bowman
Andrew Foster
Madison Stemm

With support from:

Jonathan Fortney (UCSC), Nikku Madhusudhan (Yale), Patricio Rojo (U. de Chile)

BART Development Plan

0. Fully understand Pato's code (all)
1. Split transit code into initialization and loop (Patricio)
2. Drive loop with arbitrary parameters (Patricio)
   (hope to get here by end of Pato/Patricio visit)
3. Interface transit code with MC driver (Jasmina, help from Patricio)
4. First test: H2O-only, transit, constant abundances, constant T(p)
   (good enough for NESSF!  15 Jan)
5. Upgrades required to publish (do in any order):
   many molecules' line lists (Patricio)
   eclipse geometry (Jasmina)
   CEA initial conditions/constant scaling (Jasmina)
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
   face map
   phase curve
   face map + face curve

Test after each upgrade, record output, compare output of next upgrade
with upgrade turned off to previous output to make sure we didn't
wreck anything.
