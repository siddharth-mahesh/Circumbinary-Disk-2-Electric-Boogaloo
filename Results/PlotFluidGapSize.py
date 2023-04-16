import numpy as np
import matplotlib.pyplot as plt
import sys,os

mass_ratios = np.arange(0.1,0.6,0.1)
inclinations = np.arange(0,1,0.01)*90

## edit this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","sidmahesh\OneDrive","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
data_path = os.path.join(repo_path,"Data")


# for i in range(len(mass_ratios)):
#     gapsizepath = "FluidGapSizeByInclination_q"+str(i+1)+".txt"
#     plot_label = r'$\mu$ = 0.'+str(i+1)
#     xgap = np.loadtxt(gapsizepath)
#     plt.plot(inclinations,xgap,label = plot_label)


# plt.xlabel(r'Inclination($^\circ$)')
# plt.ylabel(r'Gap Size ($a$)')
# plt.legend()
# plt.savefig('FluidGapSizeByInclination.png',dpi = 300)
# plt.show()
    
hires_mass_ratios = np.arange(0.01,0.51,0.01)
xgap = np.loadtxt("FluidGapSizeCoplanarHiRes.txt")
AL1994gap = np.loadtxt(os.path.join(data_path,"AL1994.txt"))
#AL1994mass_ratios = [0.05,0.1,0.2,0.3,0.4,0.5]

plt.plot(hires_mass_ratios,xgap,label = "Our Study")
plt.scatter(AL1994gap[:,0],AL1994gap[:,1],edgecolors='black',facecolors = 'none',label = 'AL94')

plt.xlabel(r'$\mu$')
plt.ylabel(r'$r_T/a$')
plt.legend()
plt.savefig('FluidGapSizeCoplanar.png',dpi = 300)
plt.show()
 
# eps = np.loadtxt('FluidGapSizeEpsVaryQ.txt')
# plt.plot(inclinations,np.log10(eps))
# plt.savefig('FluidEpsVaryQ.png',dpi = 300)
# plt.show()

# eps = np.loadtxt('FluidGapSizeEpsVaryE.txt')
# plt.plot(inclinations,np.log10(eps))
# plt.savefig('FluidEpsVaryE.png',dpi = 300)
# plt.show()

