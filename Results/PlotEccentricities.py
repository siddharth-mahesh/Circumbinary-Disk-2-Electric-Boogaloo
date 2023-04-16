import numpy as np
import matplotlib.pyplot as plt

for i in range(5):
    res_file_name = "EccentricityComputation"+str(i+1)+".txt"
    data = np.loadtxt(res_file_name)
    data_label = 'q = 0.'+str(i+1)
    plt.plot(data[:,0],np.log10(data[:,1]),label = data_label)

plt.xlabel("Radial Distance")
plt.ylabel("Log(Eccentricity)")
plt.legend()
plt.savefig("Eccentricities.png")
plt.show()
