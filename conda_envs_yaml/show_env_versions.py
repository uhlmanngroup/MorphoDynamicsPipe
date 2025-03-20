#python show_env_versions.py
#/nfsresearch/bwoodhams/plast_cell/2024-08-22_making_example/show_env_versions.ipynb

import os
import subprocess

directory_out = os.path.abspath('.')
directory_in = os.path.abspath('.')
print(os.listdir(directory_in))
envs = [each.replace('_dev.yml', '').replace('environment_', '') for each in os.listdir(directory_in) if each.endswith('yml')]
envs
for env in envs:
    try:
        command = 'conda list -n ' + env + ' > ' + os.path.join(directory_out, env) + '.txt'
        print(command)
        subprocess.run(command, shell = True)
    except:
        print('Failed to run command: ' + command)
        continue

#conda list -n morphody39 > /nfs/research/uhlmann/bwoodhams/plast_cell/2024-08-22_making_example/MorphoDynamicsPipe/conda_envs_yaml/morphody39.txt