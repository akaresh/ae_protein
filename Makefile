#!/bin/bash

.PHONY: install setup lint con

default: setup


install:
	conda env create -f fqa.yaml

setup:
	conda env update --file fqa.yaml --prune
	
con:
	python3 conda_env_export.py -n fqa -p python=3.8 > fqa.yaml

lint:
	pycodestyle --ignore=W1,W293,E701,E221 trainers/
	pycodestyle --ignore=W1,W293,E701,E221 fragment_generators/
	pycodestyle --ignore=W1,W292,E701,E221 visualization/