#!/bin/bash

.PHONY: install setup lint con

default: setup


install:
	conda env create -f fqa.yaml

setup:
	conda env update --file fqa.yaml --prune
	
con:
	conda env export --no-builds --from-history --name fqa > fqa.yaml

lint:
	pycodestyle --ignore=W1,W293,E701,E221 trainers/
	pycodestyle --ignore=W1,W293,E701,E221 fragment_generators/
	pycodestyle --ignore=W1,W292,E701,E221 visualization/