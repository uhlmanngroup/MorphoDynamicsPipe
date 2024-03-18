import os, sys
import numpy as np

#this_input = sys.argv[1]
#this_output = sys.argv[2]

this_input = list(snakemake.input)
this_output = snakemake.output[0]

#print("this_input is in the python file ", this_input)
#print('type ', type(this_input))

print("this_output is in the python file ", this_output)

np.save(this_output, np.zeros(1))