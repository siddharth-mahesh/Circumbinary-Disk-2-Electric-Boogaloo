import numpy as np
import matplotlib.pyplot as plt

gap_sizes = []
inclinations = np.arange(0.,1.,1./180.)*np.pi
mass_ratios = np.arange(0.1,0.6,0.1)

gap_sizes_file_name = "InclinedGapSizes.txt"
for j in range(len(inclinations)):
    inc = inclinations[j]
    gap_sizes_for_this_inc = []
    for i in range(len(mass_ratios)):
        res_file_name = "InclinedLyapExpComputation_q"+str(i+1)+"_i"+str(inc)+".txt"
        print(res_file_name)
        data = np.loadtxt(res_file_name)
        data_label = r'$\mu$ = 0.'+str(i+1)
        r = data[:,0]
        max_lyap = np.array([max(data[k,1:]) for k in range(len(data))])
        time_scales = 1/max_lyap/2/np.pi
        gap_sizes_for_this_inc.append([(i+1)*0.1,r[np.abs(np.log10(time_scales)).argmin()]])
        plt.plot(r,time_scales,label = data_label)
    plt.axhline(0,color = 'black',label = r'$P$')
    plt.yscale('log')
    plt.xlabel(r'$r/a$')
    plt.ylabel(r'$\tau/P$')
    plt.legend()
    figname = "InclinedInstabilityScalePlot_"+"i"+str(inc)+".png"
    plt.savefig(figname)
    plt.show()
    gap_sizes.append(gap_sizes_for_this_inc)


gap_sizes = np.array(gap_sizes)
#print(gap_sizes)
#np.savetxt("inclined_gap_sizes.txt",gap_sizes)

numdatainclined2pstable = np.array([[1.7991,1.7761,1.7733,1.7626,1.7373,1.7067,1.6735,1.4911,1.4714,1.4],[],[],[],[]])

print(gap_sizes)

for j in range(7):
    print("j = ",j)
    inc = (j)*15
    gaps_for_j = gap_sizes[j]
    print(gaps_for_j)
    print(gaps_for_j[:,0])
    data_label = "inc = "+str(inc)
    plt.plot(gaps_for_j[:,0],gaps_for_j[:,1],label = data_label)

plt.ylabel("Gap Size")
plt.xlabel("Mass Ratio")
plt.legend()
#plt.tight_layout()
plt.savefig("InclinedGapSizesByMassRatio.png")

plt.show()

inclinations = np.arange(0,105,15,dtype = int)

for j in range(5):
    print(gap_sizes[:,j,1])
    data_label = r"$\mu$ = 0."+str(j+1)
    plt.plot(inclinations,gap_sizes[:,j,1],label = data_label)

plt.ylabel("Gap Size")
plt.xlabel("Inclinations")
plt.legend()
#plt.tight_layout()
plt.savefig("InclinedGapSizesByInclination.png")



