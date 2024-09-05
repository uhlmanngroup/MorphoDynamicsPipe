#python show_env_versions.py
#/nfsresearch/bwoodhams/plast_cell/2024-08-22_making_example/show_env_versions.ipynb

import os
import subprocess

directory_out = os.path.abspath('.')
directory_in = os.path.abspath('../.snakemake/conda/')
print(os.listdir(directory_in))
envs = [each.replace('.yaml', '') for each in os.listdir(directory_in) if each.endswith('yaml')]
envs
for env in envs:
    command = 'conda list -p ' + os.path.join(directory_in, env) + ' > ' + os.path.join(directory_out, env) + '.txt'
    print(command)
    subprocess.run(command, shell = True)

#conda list -n morphody37 > /nfs/research/uhlmann/bwoodhams/plast_cell/2024-08-22_making_example/MorphoDynamicsPipe/conda_envs_yaml/morphody37.txt