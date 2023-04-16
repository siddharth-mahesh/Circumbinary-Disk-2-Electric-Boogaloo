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

import ResonantTorquePicture as rtp

print("initializing arrays")
## set choices for mass ratios and radio

mass_ratios = np.arange(0.1,0.5+0.1,0.1)
inclinations = np.arange(0,1,1/180)*np.pi
eccentricities = np.arange(0.0,0.5,0.01)

## test case: R = 1e4

alpha = 1e-2
chi = 1e-1

## Inclined non eccentric cases

# for i in range(len(mass_ratios)):
#     print("\n q = %.2e \n"%mass_ratios[i])
#     res_file_name = "FluidGapSizeByInclination_q"+str(i+1)+".txt"
#     res_file = os.path.join(results_path,res_file_name)
#     xgap = np.ones(len(inclinations))
#     for j in range(len(inclinations)):
#         rgap = 0
#         rgapparams = np.zeros(6)
#         incl = inclinations[j]
#         for m in range(5,0,-1):
#             rgapguess = rtp.LR_location(m,1)
#             params = [rgapguess,mass_ratios[i],incl,0.0,m,1]
#             fluidparams = [params,alpha,chi]
#             zeta = rtp.zeta_T(fluidparams)
#             if zeta > 1:
#                 if rgapguess > rgap:
#                     rgap = rgapguess
#                     rgapparams = params
                    
#         rgapguess = rtp.LR_location(2,2)
#         params = [rgapguess,mass_ratios[i],incl,0.0,2,2]
#         fluidparams = [params,alpha,chi]
#         zeta = rtp.zeta_T(fluidparams)
#         if zeta > 1:
#             if rgapguess > rgap:
#                 rgap = rgapguess
#                 rgapparams = params
#         if rgap != 0:
#             xgap[j] = rtp.find_rgap( rgapparams , [rgap+1e-10,rgap+2])
#         #xgap[j] = rgap
#     np.savetxt(res_file,xgap)

## non-inclined non eccentric case hi-res

print("non-inclined non eccentric case hi-res")
hires_mass_ratios = np.arange(0.0125,0.5125,0.0125)
xgap = np.zeros(len(hires_mass_ratios))
alpha = 1e-2
chi = 1e-1
res_file = os.path.join(results_path,"FluidGapSizeCoplanarHiResHydro.txt")

for i in range(len(hires_mass_ratios)):
    rgap = 0
    for m in range(1,6):
        rgapguess = rtp.LR_location(m,m)
        params = [rgapguess,hires_mass_ratios[i],0.0,0.0,m,m]
        fluidparams = [params,alpha,chi]
        zeta = rtp.zeta_T(fluidparams)
        if zeta > 1:
            if rgapguess > rgap:
                rgap = rgapguess
                rgapparams = params
    if rgap!=0:
        xgap[i] = rtp.find_rgap(rgapparams, [rgap+1e-10,rgap+2])
np.savetxt(res_file,xgap)

## find epsilon for each mass ratio

# print("finding epsilon for each inclination when the m = 2 transition occurs by varying the mass ratio")

# hires_mass_ratios = np.arange(0.05,0.50,0.001)
# res_file = os.path.join(results_path,"FluidGapSizeEpsVaryQ.txt")
# eps_for_inc = np.zeros(len(inclinations))
# for j in range(len(inclinations)):
#     for i in range(len(hires_mass_ratios)-1,0,-1):
#         rgapguess = rtp.LR_location(1,1)
#         params = [rgapguess,hires_mass_ratios[i],inclinations[j],0.0,1,1]
#         zeta = rtp.zeta_T(fluidparams)
#         if zeta > 1:
#             break
#     eps_for_inc[j] = 0.5 - hires_mass_ratios[i]

# np.savetxt(res_file,eps_for_inc)

# print("finding epsilon for each inclination when the m = 2 transition occurs by varying the eccentricity")

# hires_ecc = np.arange(0,0.2,0.001)
# res_file = os.path.join(results_path,"FluidGapSizeEpsVaryE.txt")
# eps_for_inc = np.zeros(len(inclinations))
# for j in range(len(inclinations)):
#     for i in range(len(hires_ecc)):
#         rgapguess = rtp.LR_location(1,1)
#         params = [rgapguess,0.5,inclinations[j],hires_ecc[i],1,1]
#         fluidparams =[params,alpha,chi]
#         zeta = rtp.zeta_T(fluidparams)
#         if zeta > 1:
#             break
#     eps_for_inc[j] = hires_ecc[i]

# np.savetxt(res_file,eps_for_inc)
    


