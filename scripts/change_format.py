import skimage
import numpy as np
import sys
this_input = sys.argv[1]
this_output = sys.argv[2]
print ('argument list', sys.argv)
#np.save('results/1.npy', np.array([1]))
np.save(sys.argv[2], np.array([1]))
#myarray = 

#def do_something(data_path, out_path, threads, myparam):
    # python code
#    return out_path

#do_something(snakemake.input[0], snakemake.output[0], snakemake.threads, snakemake.config["myparam"])
print(1)