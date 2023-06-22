import numpy as np
import matplotlib.pyplot as plt
import sys,os

## edit this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","sidmahesh/OneDrive","Documents\GitHub\Circumbinary-Disk-2-Electric-Boogaloo")

## no need to edit these paths
data_path = os.path.join(repo_path,"Data")

# Define MP Mass Ratios and Inclinations
mass_ratios = ["0.0099","0.099","0.2","0.3","0.4","0.5"]
mr_labels = ["1100","110","14","37","23","11"]
alpha_labels = ["0.010","0.030","0.003"]
alpha_labels_out = ["0010","0030","0003"]
alphas = [r'$1\times10^{-2}$',r'$3\times10^{-2}$',r'$3\times10^{-3}$']
mass_ratios_hires = np.arange(0.0125,0.5125,0.0125)

for j in range(len(alpha_labels)):
    gapfile_name = os.path.join(data_path,"R2P_alpha_"+alpha_labels_out[j]+".txt")
    gapfile = open(gapfile_name,"w")
    for i in range(len(mass_ratios)):
        datfile = np.loadtxt(os.path.join(data_path,"Paper1_Viscous_Data",mr_labels[i]+"_"+alpha_labels[j]+"_8.dat"))
        radii = datfile[:,3]
        cutoff = np.argmin(np.abs(radii - 3))
        densities = datfile[:,4]
        maxdense = np.max(densities)
        sig2p = 0.02*maxdense
        #if mr_labels[i] =="5":
        #    if inc_labels_125[j] == "180":
        #        check = 1
        #        #print(np.abs(densities[2:cutoff] - .02*maxdense))
        sig2p_index = np.argmin(np.abs(densities[2:cutoff] - .02*maxdense)) + 2
        if densities[sig2p_index] > sig2p:
            sigma_plus, r_plus = densities[sig2p_index] , radii[sig2p_index]
            sigma_minus, r_minus = densities[sig2p_index - 1] , radii[sig2p_index - 1]
            
        else:
            sigma_plus, r_plus = densities[sig2p_index+1] , radii[sig2p_index+1]
            sigma_minus, r_minus = densities[sig2p_index] , radii[sig2p_index]
            
        r2p = r_minus + (r_plus - r_minus)*(sig2p - sigma_minus)/(sigma_plus - sigma_minus)
        gapfile.write(mass_ratios[i] + "\t"  + str(r2p) + "\t"+ str(maxdense) + "\t"+str(sigma_minus)+"\t"+str(sigma_plus)+"\n")
    gapfile.close()
    r2pdata = np.loadtxt(gapfile_name)
    plt.scatter(r2pdata[:,0],r2pdata[:,1],label = alphas[j])
    plt.scatter(r2pdata[:,0],r2pdata[:,1],label = alphas[j])
    if (alpha_labels_out[j] == "0010"):
        rl = np.loadtxt(os.path.join(repo_path,"Results","coplanar_gap_sizes.txt"))
        rt = np.loadtxt(os.path.join(repo_path,"Results","FluidGapSizeCoplanarHiResHydro.txt"))
        plt.plot(mass_ratios_hires,rl*r2pdata[-1,1]/rl[-1],label = r"$\bar{r}_L$")
        plt.plot(mass_ratios_hires,rt*r2pdata[-1,1]/rt[-1],label = r"$\bar{r}_T$")

plt.legend()
plt.title("8000BR")
plt.savefig(os.path.join(repo_path,"Results","GapSizeCoplanarCompareViscosities_8000.png"),dpi = 300)
plt.show()
    
    