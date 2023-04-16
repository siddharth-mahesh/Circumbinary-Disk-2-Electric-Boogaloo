print("imports")
import sys,os
import numpy as np

print("setting up paths")
## edit only this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","siddh","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
formulae_path = os.path.join(repo_path,"Formulae")
sys.path.insert(0,formulae_path)
data_path = os.path.join(repo_path,"Data\RP1981.txt")
results_path = os.path.join(repo_path,"Results")

print("reading in data")
## read in the RP1981 data

rp1981_data = np.loadtxt(data_path)
mass_ratios = rp1981_data[:,0]
jacobi_constants = rp1981_data[:,4]
avg_radii = 0.5*(rp1981_data[:,1] + rp1981_data[:,2])

print("primary computation")
## compute the jacobi constants

import JacobiConstant as cj

res_file = open(os.path.join(results_path,"JacobiConstantComparisons.txt"),"w")
res_index_file = open(os.path.join(results_path,"JacobiConstantComparisonsLabels.txt"),"w")
res_index_file.write("mass_ratio \t r_avg \t jacobi_constant_comp \t jacobi_constant_rp1981 \t relative_error")
res_index_file.close()
mmax = 3

for i in range(len(mass_ratios)):
    cj_params = [avg_radii[i],mass_ratios[i],mmax]
    cj_comp = cj.pert_C_j(cj_params)
    rel_err = np.abs(cj_comp/jacobi_constants[i] - 1)
    res_file.write("%f \t %f \t %f \t %f \t %f \n"%(mass_ratios[i],avg_radii[i],cj_comp,jacobi_constants[i],rel_err))

print("done")
res_file.close()
