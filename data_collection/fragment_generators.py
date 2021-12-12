#!/usr/bin/env python

import sys

from molmass import Formula
import pandas as pd

alphabet = [
	"G",
	"P",
	"D",
	"E",
	"K",
	"R",
	"H",
	"S",
	"T",
	"N",
	"Q",
	"A",
	"M",
	"Y",
	"W",
	"V",
	"I",
	"L",
	"F",
	"C"
]


aa_dict = {
	"ala" : "A", 
	"arg" : "R",
	"asn" : "N",
	"asp" : "D",
	"cys" : "C",
	"gln" : "Q",
	"glu" : "E",
	"gly" : "G",
	"his" : "H",
	"ile" : "I",
	"leu" : "L",
	"lys" : "K",
	"met" : "M",
	"phe" : "F",
	"pro" : "P",
	"ser" : "S",
	"thr" : "T",
	"trp" : "W",
	"tyr" : "Y",
	"val" : "V"
}


aa_sc_atoms = {
	"A" : ["CB"], 
	"R" : ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
	"N" : ["CB", "CG", "OD1", "ND2"],
	"D" : ["CB", "CG", "OD1", "OD2"],
	"C" : ["CB", "SG"],
	"Q" : ["CB", "CG", "CD", "OE1", "NE2"],
	"E" : ["CB", "CG", "CD", "OE1", "OE2"],
	"G" : ["CA"],
	"H" : ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
	"I" : ["CB", "CG1", "CG2", "CD1"],
	"L" : ["CB", "CG", "CD1", "CD2"],
	"K" : ["CB", "CG", "CD", "CE", "NZ"],
	"M" : ["CB", "CG", "SD", "CE"],
	"F" : ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
	"P" : ["CB", "CG", "CD", "N", "CA"],
	"S" : ["CB", "OG"],
	"T" : ["CB", "OG1", "CG2"],
	"W" : ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
	"Y" : ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
	"V" : ["CB", "CG1", "CG2"]
}


def ca_fragment(struct, chain, ids, dssp):
	pos = []
	
	for i in ids:
		if i not in struct[0][chain]: return None
		if "CA" not in struct[0][chain][i]: return None
		
		residue = struct[0][chain][i]
		if residue.get_resname().lower() not in aa_dict: return None
		
		grp = residue["CA"].get_coord().tolist()
		pos.append(grp)
	
	return pos


def bb_fragment(struct, chain, ids, dssp):
	bbatms = ['N', 'CA', 'C', 'O']
	pos = []
	
	for i in ids:
		if i not in struct[0][chain]: return None
		
		residue = struct[0][chain][i]
		if residue.get_resname().lower() not in aa_dict: return None
		
		for bb in bbatms:
			if bb not in struct[0][chain][i]: return None
			
			grp = residue[bb].get_coord().tolist()
			pos.append(grp)
	
	return pos


def res_cen(dict_positions):
	#print(dict_positions)
	total_mass = 0
	mx = 0
	my = 0
	mz = 0
	for keys, values in dict_positions.items():
		total_mass += (Formula(values[-1]).mass)
		mx += ((Formula(values[-1]).mass)*values[0])
		my += ((Formula(values[-1]).mass)*values[1])
		mz += ((Formula(values[-1]).mass)*values[2])
	return [mx/total_mass, my/total_mass, mz/total_mass]


def bbcen_fragment(struct, chain, ids, dssp):
	bbatms = ["N", "CA", "C", "O"]
	pos = []
	
	for i in ids:
		if i not in struct[0][chain]: return None
		
		residue = struct[0][chain][i]
		if residue.get_resname().lower() not in aa_dict: return None
		
		for bb in bbatms:
			if bb not in struct[0][chain][i]: return None
			
			grp = residue[bb].get_coord().tolist()
			pos.append(grp)
		
		resname = aa_dict[residue.get_resname().lower()]
		if resname == "G":
			pos.append(struct[0][chain][i]["CA"].get_coord().tolist())
			continue
		
		sc_atms_check = {k:0 for k in aa_sc_atoms[resname]}
		side_coords = dict()
		for atm in residue.get_atoms():
			if resname != "P" and resname != "G":
				if atm.get_id() in bbatms: continue
			
			if atm.get_id()[0] == 'H': continue
			if atm.get_id()[0] == 'D': continue
			
			if atm.get_id() not in sc_atms_check: continue
			assert(sc_atms_check[atm.get_id()] == 0)
			assert(atm.get_id() not in side_coords)
			
			side_coords[atm.get_id()] = [
				atm.get_coord()[0],
				atm.get_coord()[1],
				atm.get_coord()[2],
				atm.element]
			sc_atms_check[atm.get_id()] += 1
		
		for k,v in sc_atms_check.items():
			if v != 1: return None
		
		ressc_centroid = res_cen(side_coords)
		pos.append(ressc_centroid)
	return pos


def bbsc_fragment(struct, chain, ids, dssp):
	pass


























