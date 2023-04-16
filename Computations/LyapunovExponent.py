print("imports")
import sys,os
import numpy as np

print("setting up paths")
## edit only this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","sidmahesh\OneDrive","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
formulae_path = os.path.join(repo_path,"Formulae")
sys.path.insert(0,formulae_path)
results_path = os.path.join(repo_path,"Results")

print("initializing arrays")
## set choices for mass ratios and radio

mass_ratios = np.arange(0.0,0.5+0.1,0.1)
inclinations = np.arange(0,1,1/180)*np.pi
eccentricity = np.arange(0.0,0.5,0.1)
avg_radii = np.arange(1.10,3,0.001)

print("primary computation")
## compute the lyapunov exponents

import RotatingStabilityMatrix as rsm
import InclinedStabilityMatrix as ism
from scipy.linalg import eigvals

## params get parsed to other functions and potentials as
## r, q, i, ecc, mmin, Nmin, mmax, Nmax


mmax = 3
mmin = 1
nmax = 5 
nmin = 1

res_index_file = open(os.path.join(results_path,"LyapExpComputationLabels.txt"),"w")
res_index_file.write("r_avg \t l0 \t l1 \t l2 \t l3")
res_index_file.close()

test_file_name = "test_keplerian.txt"
test_file = open(os.path.join(results_path,test_file_name),"w")

inclined = 0
eccentric = 0

coplanar_mass_ratios = np.arange(0.0125,0.5125,0.0125)
if inclined == 0:
    for i in range(len(coplanar_mass_ratios)):
        q = coplanar_mass_ratios[i]
        res_file_name = "CoplanarLyapExpComputation_q_"+str(int(1000*q))+".txt"
        res_file = open(os.path.join(results_path,res_file_name),"w")
        print("mass ratio : %f"%(q))
        for j in range(len(avg_radii)):
            #print("r = ", avg_radii[j])
            params = [avg_radii[j],coplanar_mass_ratios[i],0.0,0.0,mmin,nmin,mmax,nmax]
            stability_mat = rsm.K(params)
            #print(stability_mat)
            lyapexps = eigvals(stability_mat)
            #print(lyapexps)
            if i == 0:
                test_file.write("%f \t %f \t %f \t %f \t %f \n"%(avg_radii[j],lyapexps[0].imag,lyapexps[1].imag,lyapexps[2].imag,lyapexps[3].imag))
            res_file.write("%f \t %f \t %f \t %f \t %f \n"%(avg_radii[j],lyapexps[0].real,lyapexps[1].real,lyapexps[2].real,lyapexps[3].real))
        res_file.close()

# if inclined == 1:
#     for i in range(len(mass_ratios)):
#         q = mass_ratios[i]
#         print("mass ratio : %f"%(q))
#         for j in range(len(inclinations)):
#             inc = inclinations[j]
#             inc_in_deg = j
#             print("inclination : %f"%(inc_in_deg))
#             res_file_name = "InclinedLyapExpComputation_q"+str(int(10*q))+"_i"+str(inc_in_deg)+".txt"
#             res_file = open(os.path.join(results_path,res_file_name),"w")
#             for k in range(len(avg_radii)):
#                 #print("r = ", avg_radii[k])
#                 params = [avg_radii[k],mass_ratios[i],inclinations[j],0.0,mmin,nmin,mmax,nmax]
#                 stability_mat = ism.K(params)
#                 #print(stability_mat)
#                 lyapexps = eigvals(stability_mat)
#                 #print(lyapexps)
#                 if i == 0 and j == 0:
#                     test_file.write("%f \t %f \t %f \t %f \t %f \n"%(avg_radii[k],lyapexps[0].imag,lyapexps[1].imag,lyapexps[2].imag,lyapexps[3].imag))
#                 res_file.write("%f \t %f \t %f \t %f \t %f \n"%(avg_radii[k],lyapexps[0].real,lyapexps[1].real,lyapexps[2].real,lyapexps[3].real))
#             res_file.close()
#             if i == 0:
#                 break

test_file.close()
print("done")