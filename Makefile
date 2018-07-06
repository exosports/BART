# Makefile for BART
#
# `make` - Compile submodules and manuals
#
# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

all: make_man make_transit make_mccubed
clean: clean_userman clean_transit clean_mccubed

make_man:
	@echo "\nCompiling user manual..."
	@cd doc/BART_user_manual/ && pdflatex BART_user_manual.tex
	@echo "\nFinished compiling user manual."

make_transit:
	@echo "\nBuilding transit module..."
	@cd modules/transit/ && make
	@echo "\nFinished building transit module."

make_mccubed:
	@echo "\nBuilding MCcubed module..."
	@cd modules/MCcubed && make
	@echo "\nFinished building MCcubed module."

clean_userman:
	@echo "\nRemoving user manual..."
	@cd doc/BART_user_manual/ && rm -f BART_user_manual.pdf
	@echo "\nDone."

clean_transit:
	@echo "\nRemoving transit module binaries..."
	@cd modules/transit/ && make clean
	@echo "\nDone."

clean_mccubed:
	@echo "\nRemoving MCcubed module binaries..."
	@cd modules/MCcubed/ && make clean
	@echo "\nDone."
