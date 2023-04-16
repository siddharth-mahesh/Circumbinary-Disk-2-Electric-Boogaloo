import numpy as np
import matplotlib.pyplot as plt

res_file_name = "test_keplerian.txt"
data = np.loadtxt(res_file_name)
r = data[:,0]
lyap_comp = data[:,2]
lyap_true = 1/(r**1.5)
e_rel = np.abs(lyap_comp/lyap_true - 1)
plt.plot(r,np.log10(e_rel),label = "relative error")


plt.xlabel("Radial Distance")
plt.ylabel("Log(Relative Error)")
plt.legend()
plt.savefig("TestKeplerian.png")
plt.show()

