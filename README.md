BART
====

Bayesian Atmospheric Radiative Transfer fitting code

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

Apparently GIThub's Markdown doesn't support tables yet...

|------------------|--------|-------|--------|----|----------------|  
|Component		    |Code	 |Doc	   |Package	|Git |Centralgit(here)|  
|------------------|--------|-------|--------|----|----------------|  
|comm loop		    |dev	    |-      |-	      |y	  |-               |  
|DEMC 		       |done	 |-      |-	      |-	  |-               |  
|TEA			       |done    |dev	   |-	      |y	  |-               |  
|input converter   |ready   |incode	|-	      |-	  |-               |  
|transit		       |dev	    |thesis	|-	      |y2  |outdated        |  
|lineread	       |dev	    |-      |-	      |y	  |outdated        |  
|output converter  |ready	 |incode	|-	      |-	  |-               |  
|------------------|--------|-------|--------|----|----------------|  
