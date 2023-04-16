import numpy as np
import sys,os
import matplotlib.pyplot as plt

## edit this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","sidmahesh\Onedrive","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
data_path = os.path.join(repo_path,"Data")

gap_sizes = []
linestyles = [(0,(5,1)),(5,(10,3)),(0,(5,5)),'dashed','dashdot']

for i in range(5):
    res_file_name = "InclinedLyapExpComputation_q"+str(i+1)+"_i0.txt"
    data = np.loadtxt(res_file_name)
    data_label = r'$\mu$ = 0.'+str(i+1)
    r = data[:,0]
    max_lyap = np.array([max(data[j,1:]) for j in range(len(data))])
    #time_scales = np.log10(1/max_lyap/2/np.pi)
    time_scales = 1/max_lyap/2/np.pi
    #gap_sizes.append([(i+1)*0.1,r[np.abs(time_scales).argmin()]])
    plt.plot(r,time_scales,label = data_label,linestyle = linestyles[i])


#plt.axhline(0,color = 'black',label = r'$\delta = 1$')
plt.axhline(1,color = 'black',label = r'$P$')
plt.yscale('log')
plt.ylim(1e-2,1e1)
plt.xlabel("$r/a$")
plt.ylabel(r'$\tau/P$')
plt.legend()
plt.savefig("InstabilityScale.png",dpi = 300)
plt.show()


coplanar_mass_ratios = np.arange(0.0125,0.5125,0.0125)
gap_sizes = np.zeros(len(coplanar_mass_ratios))
for i in range(len(coplanar_mass_ratios)):
    res_file_name = "CoplanarLyapExpComputation_q_"+str(int(1000*coplanar_mass_ratios[i]))+".txt"
    data = np.loadtxt(res_file_name)
    r = data[:,0]
    max_lyap = np.array([max(data[j,1:]) for j in range(len(data))])
    time_scales = np.log10(1/max_lyap/2/np.pi)
    gap_sizes[i] = r[np.abs(time_scales).argmin()]
    
np.savetxt("coplanar_gap_sizes.txt",gap_sizes)
# plt.plot(coplanar_mass_ratios,gap_sizes)
# plt.ylabel(r'$r_\mathrm{L}$')
# plt.xlabel(r'$\mu$')
# plt.savefig("GapSizesRL.png")
# plt.show()


numdatacoplanar2p = np.array([1.799,1.782,1.744,1.676,1.500,1.188])
numdatacoplanardT = np.array([2.078,2.070,2.048,1.992,1.861,1.402])
numdataq = ([.5,.4,.3,.2,.09,.0099])
hires_mass_ratios = np.arange(0.0125,0.5125,0.0125)
r_T = np.loadtxt("FluidGapSizeCoplanarHiResHydro.txt")

plt.plot(coplanar_mass_ratios,gap_sizes,color = 'black', label = r'$r_\mathrm{L}$',linestyle = 'dashed')
plt.plot(hires_mass_ratios,r_T*gap_sizes[-1]/r_T[-1], color = 'red',label = r'$\bar{r}_\mathrm{T}$',linestyle = (0,(5,5)))
plt.scatter(numdataq,gap_sizes[-1]*numdatacoplanar2p/numdatacoplanar2p[0],edgecolors = 'black',facecolors = 'none',label = r'$\bar{r}_\mathrm{2\%}$',marker = 'o')
plt.scatter(numdataq,gap_sizes[-1]*numdatacoplanardT/numdatacoplanardT[0],edgecolors = 'red',facecolors = 'none',label = r'$\bar{r}_\mathrm{dT}$',marker = '^')
plt.ylabel(r'$r/a$')
plt.xlabel(r'$\mu$')
plt.legend()
plt.savefig("GapSizesCoplanarRLR1p.png",dpi=300)
plt.show()


