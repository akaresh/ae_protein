#!/usr/bin/env python

import argparse
import os
import sys

parser = argparse.ArgumentParser(description='creating turn labels')
parser.add_argument('--pdbs', required=True, type=str,
	metavar='<path>', help='path to chains directory')
parser.add_argument('--bt18', required=True, type=str,
	metavar='<path>', help='path to BetaTurns18.py2')
parser.add_argument('--dssp', required=True, type=str,
	metavar='<path>', help='path to binary for dssp')

arg = parser.parse_args()

for pdb in os.listdir(arg.pdbs):
	if os.path.isfile(os.path.join(arg.pdbs, pdb)):
		if pdb.endswith('.dssp'):
			splits = pdb.split('.')
			pdbid = splits[0]
			#print(pdbid)
			input_pdb = os.path.join(arg.pdbs, pdbid+'.pdb')
			#print(input_pdb)
			assert(os.path.isfile(input_pdb))
			output_bt18 = os.path.join(arg.pdbs, pdbid+'.bt18.out')
			#print(output_bt18)
			cmd = "python {bt18} -i {input} --dssp {dssp} -noG --quiet\
			--no-seq-format --no-comments -o {output}".format(
				bt18=arg.bt18,
				dssp=arg.dssp,
				input=input_pdb,
				output=output_bt18)
			print(cmd)
			os.system(cmd)
			#sys.exit()
"""
Options
1. Our project clusters beta turns, 4 residue turns
	- identified using BetaTurnTool
	- Pros: we don't need to worry with making turns, just get from 
	BetaTurnTool output
	- Cons: The focus on defining beta turns as any 4 residue segment where
	middle two residues are not H,E and i->i+3 CA distance < 7 makes overlapping
	4 residue "turn" segments. They are only part of the turn. 

1. just use the collected beta turns from a 50% filtering of the PDB
	- cons: overlapping 4 residue turns
	- cons: we adding context, will we get any clusters in the latent space
	- with our approach, the clusters/latent space that we will learn its
	going to be able to reconstruct the fragment
	- well defined clusters, but we will learn a network that has an informative
	latent space, informative enough to reconstruct
	-
just focus on beta turns --> more turns, ~30 turns per protein
other types of turns --> 
the way we currently do it
we just DSSP T labels, no other labels can comprise a turn
we wouldnt use DSSP to label beta turns
2,3,4,5,6
no
2 tracks
1. beta turn definitions from dunbrack
2. DSSP only, only using residues with the T label
3. DSSP only for beta turns --> yes, it should be nearly the same as Dunbrack
4. Dunbrack identities the beta turns, and uses their own labeling/new labeling to 
label the turns, 14 new types, 4 residues, middle 2/3 are not H/E, and 1->4 is < 7A


so to use context or not
latent space that all the turns live
populate that space with the turns, their image in that space
those points will get clustered by some other clustering method, 
PCA or DBSCAN --> then you get the clusters and cluster types
if we add context, there will be many points that cant get clusters
but, but network knows how to reconstuct them
the network has picked on structural features that makes turns distinct

in dbscan for example, there is mechanism for points to not be clustered
PCA is used in clustering, but its not a clustering algorithm --> lower dimension projection
dbscan is clustering 
latent space -> pca -> dbscan
latent space -> dbscan

matching up, focusing on beta turns could be next step or a test of the dunbrack method
context or no context, and we take distance matrix, and learn on that
first, can we get a low reconstruction loss
if so, what does latent space look like
and does map/match with the densities/types from Dunbrack 18?

4 fragment types
3 levels of context
pca first or not
(hyper parameter selection of CNN AE)

post processing 
we dont allow any turns that overlap
if any two turns starting index difference is 1, then you throw both of them out

thursday just focus on base model just applying pca to the coordiantes
dbscan








 
"""
