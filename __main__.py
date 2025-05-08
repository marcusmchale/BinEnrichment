#!/usr/bin/env python
import argparse
from src import InputHandler
from src import TestHandler


parser = argparse.ArgumentParser()
parser.add_argument('-m','--mapping', required=True)
parser.add_argument('-o','--output', required=True)
parser.add_argument('-a','--alpha', default=0.05, type=float)
# Files should be simple with each line containing a value that corresponds to an identifier in the mapping file
parser.add_argument('-bf', '--background-file')
parser.add_argument('-tf', '--target-file')
parser.add_argument('-uf', '--up-file')
parser.add_argument('-df', '--down-file')
# Lists should be bash arrays of identifiers found in the mapping file
parser.add_argument('-bl', '--background-list', nargs='+')
parser.add_argument('-tl', '--target-list', nargs='+')
parser.add_argument('-ul', '--up-list', nargs='+')
parser.add_argument('-dl', '--down-list', nargs='+')
# If both file and list are provided, they will be merged

args = vars(parser.parse_args())
parsed_input = InputHandler(args)
TestHandler(
	parsed_input.mapping,
	parsed_input.output,
	parsed_input.background,
	diff=parsed_input.target,
	up=parsed_input.up,
	down=parsed_input.down,
	alpha=parsed_input.alpha
).perform_test()
