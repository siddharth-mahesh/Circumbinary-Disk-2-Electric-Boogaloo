print("imports")
import sys,os
import numpy as np

print("setting up paths")
## edit only this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","siddh","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
formulae_path = os.path.join(repo_path,"Formulae")
sys.path.insert(0,formulae_path)
results_path = os.path.join(repo_path,"Results")

print("initializing arrays")
## set choices for mass ratios and radio

mass_ratios = np.arange(0.1,0.5+0.1,0.1)
avg_radii = np.arange(1.3,3,0.05)

print("primary computation")
## compute the jacobi constants

import Eccentricities as ecc
mmax = 3
res_index_file = open(os.path.join(results_path,"EccentricityComputationLabels.txt"),"w")
res_index_file.write("r_avg \t eccentricity")
res_index_file.close()

for i in range(len(mass_ratios)):
    q = mass_ratios[i]
    res_file_name = "EccentricityComputation"+str(int(10*q))+".txt"
    res_file = open(os.path.join(results_path,res_file_name),"w")
    print("mass ratio : %f"%(q))
    for j in range(len(avg_radii)):
        ecc_params = [avg_radii[j],mass_ratios[i],mmax]
        ecc_comp = ecc.avg_eccentricity(ecc_params)
        res_file.write("%f \t %f \n"%(avg_radii[j],ecc_comp))
    res_file.close()

print("done")
