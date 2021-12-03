#!/usr/bin/python3 

import argparse
import csv
import json
import os
import sys

def make_ssindex(dssp):
	ssindex = dict()
	with open(dssp, 'r') as fp:
		lines = fp.readlines()
		
		begin = False
		for line in lines:
			if '#  RESIDUE AA' in line:
				begin = True
				continue
			if begin:
				if line[11] == ' ': continue
				
				line = line.rstrip()
				row = line.split()
				
				res = row[1]
				aa  = line[13]
				ss  = line[16]
				if ss == ' ': ss = 'C'
				assert(res not in ssindex)
				ssindex[res] = (aa, ss)
	return ssindex

def turn_check(start, end, seq, ss, ssind):
	inds = list(range(start, end+1))
	for i, aa, si in zip(inds, seq, ss):
		if str(i) not in ssind: return False
		if ssind[str(i)][0] != aa: return False
		if ssind[str(i)][1] != si: return False
	
	return True

parser = argparse.ArgumentParser(description='creating bt18 processed files')
parser.add_argument('--pdbs', required=True, type=str,
	metavar='<path>', help='path to chains directory')

arg = parser.parse_args()

for pdb in os.listdir(arg.pdbs):
	if os.path.isfile(os.path.join(arg.pdbs, pdb)):
		if pdb.endswith('.bt18.out'):
			print(pdb)
			
			splits = pdb.split('.')
			pdbid = splits[0]
			dssp_file = pdbid+'.dssp'
			ssdic = make_ssindex(os.path.join(arg.pdbs, dssp_file))
			
			file = os.path.join(arg.pdbs, pdb)
			
			turns = {}
			groups = {}
			firsts = []
			turn_id = 0
			group_code = 0
			with open(file, 'rt', newline='') as fp:
				header = fp.readline()
				cols = header.split()
				
				for row in fp.readlines():
					data = row.split()
					dic = {k: v for k, v in zip(cols, data)}
					try:
						first = int(dic['res1'])
						four  = int(dic['res4'])
					except:
						continue
					
					if four - first != 3: continue
					
					if turn_check(
						int(dic['res1']),
						int(dic['res4']),
						dic['aa1234'],
						dic['ss1234'],
						ssdic):
						
						turn_id += 1
						turns[turn_id] = {
							'res1'   : int(dic['res1']),
							'res4'   : int(dic['res4']),
							'seq'    : dic['aa1234'],
							'ss'     : dic['ss1234'],
							'type1'  : dic['type1'],
							'symbol' : dic['letter1']
						}
					else: continue
					
					if len(firsts) == 0:
						firsts.append((int(dic['res1']), turn_id))
						continue
					else:	
						for t1 in firsts:
							if abs(int(dic['res1']) - t1[0]) == 1:
								if t1[1] in groups: groups[turn_id] = groups[t1[1]]
								else:
									groups[t1[1]] = group_code + 1
									groups[turn_id] = groups[t1[1]]
			
			for grps in groups:
				del turns[grps]
				
			if turns:
				turns_out = os.path.join(arg.pdbs, pdbid+'.turns.json')
				dssp_out = os.path.join(arg.pdbs, pdbid+'.dssp.json')
				with open(turns_out, 'w') as fp:
					json.dump(turns, fp)
				
				with open(dssp_out, 'w') as fp:
					json.dump(ssdic, fp)
			
"""

write out the turns .json

"""			
			
			