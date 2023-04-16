import numpy as np
import matplotlib.pyplot as plt
import sys,os

## edit this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","sidmahesh/OneDrive","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
data_path = os.path.join(repo_path,"Data")

## the dataset of azimuthal modes are labelled
## ab_yy.dat -> yy = inclination in deg 2 digits, ab = mass ratio q a:b


mass_ratios = [0.1,0.2,0.3,0.4,0.5]
inclinations = [0,15,30,45,60,75,90]
abvals = ['10','14','23','37','11']
yyvals = ['00','15','30','45','60','75','90']

for i in range(len(abvals)):
    for j in range(len(yyvals)):
        file_name = os.path.join(data_path,abvals[i]+'_'+yyvals[j]+'.dat')
        data = np.loadtxt(file_name)
        max_density_index = np.argmax(data[:,1])
        r_max = data[max_density_index,0]
        r_limit = r_max+1
        max_index = np.argmin(np.abs(data[:,0] - r_limit))
        C_1 = np.sqrt(data[10:max_index,2]**2 + data[10:max_index,3]**2)/data[10:max_index,1]
        C_2 = np.sqrt(data[10:max_index,4]**2 + data[10:max_index,5]**2)/data[10:max_index,1]
        C_sum_modes = np.zeros(len(data[10:max_index]))
        for mode in range(2,8):
            C_sum_modes += np.sqrt(data[10:max_index,2*mode]**2 + data[10:max_index,2*mode+1]**2)/data[10:max_index,1]
        dataset_name = 'C_m_'+abvals[i]+'_'+yyvals[j]
        plt.plot(data[10:max_index,0],np.log10(C_1),label = 'm=1')
        plt.plot(data[10:max_index,0],np.log10(C_2),label = 'm=2')
        plt.plot(data[10:max_index,0],np.log10(C_sum_modes), label = r'sum m $\geq$ 2')
        plt.axvline(r_max)
        plt.ylim(-5,0)
        plt.xlabel('r')
        plt.ylabel(r'$C_m$')
        plt.legend()
        plt.savefig(dataset_name+'.png',dpi = 300)
        plt.show()
        C_m_vals = np.zeros([len(C_1),4])
        C_m_vals[:,0] = data[10:max_index,0]
        C_m_vals[:,1] = C_1
        C_m_vals[:,2] = C_2
        C_m_vals[:,3] = C_sum_modes
        np.savetxt(dataset_name+'.txt',C_m_vals)
        
