
import argparse

from data_collection_library import make_CAframe, make_bbframe, make_bbcen, make_bbsc


parser = argparse.ArgumentParser(description='CNN model')
parser.add_argument('--file1', required=True, type=str,
	metavar='<path>', help='path to fasta file1')
parser.add_argument('--file0', required=True, type=str,
	metavar='<path>', help='path to fasta file0')
parser.add_argument('--split', required=False, type=float, default = 0.2,
	metavar = '<float>', help = 'split size')
parser.add_argument('--epoch', required=False, type=int, default = 10,
	metavar='<int>', help='number of epochs')
parser.add_argument('--batch', required=False, type=int, default = 2,
	metavar='<int>', help='batch size')
parser.add_argument('--seed', required=False, type=int,
	metavar='<int>', help='random seed')
arg = parser.parse_args()


def make_fragment_frame(master, size, t):
	if t == 'CA':
		return make_CAframe(master, size)
	elif t == 'bb':
		return make_bbframe(master, size)
	elif t == 'bb+cen':
		return make_bbcen(master, size) # separate function to calculate a centroid
	elif t == 'bb+sc':
		return make_bbsc(master, size) # function to calculate cb-extended for glycine
	else:
		raise Exception(f'un-suppported fragment tyoe: {t}')
		return None
