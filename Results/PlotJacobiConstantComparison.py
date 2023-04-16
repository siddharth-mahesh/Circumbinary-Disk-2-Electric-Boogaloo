import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("JacobiConstantComparisons.txt")
mass_ratio = data[:,0]
jacobi_constant_comp = data[:,2]
jacobi_constant_rp1981 = data[:,3]
relative_error = data[:,4]

plt.plot(mass_ratio,jacobi_constant_comp,label = 'epicyclic')
plt.plot(mass_ratio,jacobi_constant_rp1981,label = 'RP1981')
plt.xlabel(r"Mass ratio: $m_2/M$")
plt.ylabel(r"Jacobi constant: $C_j$ (Geomtrized) ")
plt.legend()
plt.savefig("JacobiConstantComparison.png")
plt.show()

plt.plot(mass_ratio,np.log10(relative_error))
plt.xlabel(r"Mass ratio: $m_2/M$")
plt.ylabel(r"Log Relative Error in Jacobi Constant")
plt.savefig("JacobiConstantError.png")
plt.show()
