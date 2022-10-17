from os import path
import sys, json, time, os

# get directories
BASE_DIR = path.dirname(path.dirname(path.abspath(__file__)))
APPLY_SRC = os.path.join('src', 'model_apply.jl')

print('\nComputing sequence...')

# call julia from python
os.system('julia --project=./src/RouteSequencing -e \'include("' + APPLY_SRC + '"); run_all("' + BASE_DIR + '")\'')

print('\nDone!')
