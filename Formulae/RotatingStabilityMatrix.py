import Potentials as pt
from scipy.linalg import eig
import numpy as np

## define the stability matrix K for coplanar case
## params get input into this module in function K as
## r, q, i, ecc, mmin, Nmin, mmax, Nmax
## params get parsed to other functions and potentials as
## r, q, i, ecc, m, N

def unpert_sol(params):
    r0 = params[0]
    return np.array([r0,r0**(-1.5)-1,0,np.sqrt(r0)])

def pert_sol(params):
    m = params[4]
    phi = pt.compute_eccentric_potential(params)
    dphi = pt.compute_eccentric_potential_d1(params)
    backg_sol = unpert_sol(params)
    r0 , omega0 = backg_sol[0] , backg_sol[1]
    D = (omega0 + 1)**2 - (m*omega0)**2
    return [-(2.*phi*(1 + 1/omega0)/r0 + dphi)/np.abs(D),-phi/omega0]

def K0(sol_backg):
    r0 , l0 = sol_backg[0] , sol_backg[3]
    r0inv = 1/r0
    r0_m2 = r0inv*r0inv
    r0_m3 = r0_m2*r0inv
    mat_0 = np.array([[0.,0.,1.,0.],[-2.*l0*r0_m3,0.,0.,r0_m2],[-r0_m3,0.,0.,2*l0*r0_m3],[0.,0.,0.,0.]])
    #print(mat_0)
    return mat_0

## define the modewise perturbed solution matrix

def K1(sol_pert,sol_backg,params):
    m = params[4]
    #print("m = ", m)
    r0 , l0 = sol_backg[0] , sol_backg[3]
    r0_inv = 1/r0
    r0_m2 = r0_inv*r0_inv
    r0_m3 = r0_m2*r0_inv
    r0_m4 = r0_m2*r0_m2
    phi = pt.compute_eccentric_potential(params)
    dphi = pt.compute_eccentric_potential_d1(params)
    ddphi = pt.compute_eccentric_potential_d2(params)
    #print(phi,'\n',dphi,'\n',ddphi)
    r1 , l1 = sol_pert[0],sol_pert[1]
    #print("r1 = ", r1/np.sqrt(2))
    #print("l1 = ", l1/np.sqrt(2))
    #print(r1,'\n',l1)
    Kpert = np.zeros([4,4])
    Kpert[1][0] = -2.*(l1*r0_m3 - 3.*l0*r1*r0_m4)
    #print("K_p_1,0 = " , Kpert[1][0]/np.sqrt(2))
    Kpert[1][3] = -2.*r1*r0_m3
    Kpert[2][0] = -6.*l0*l1*r0_m4 + 6.*r1*r0_m4 - ddphi
    Kpert[2][1] = m*dphi
    Kpert[2][3] = 2.*(l1*r0_m3 - 3.*l0*r1*r0_m4)
    Kpert[3][0] = m*dphi
    Kpert[3][1] = m*m*phi
    #print(Kpert)
    return Kpert
## define the averaged stability matrix

def K(params):
    mmin , mmax = params[4] , params[6]
    sol_backg = unpert_sol(params)
    mat = K0(sol_backg)
    #print("K0 = ", mat)
    for m in range(mmin,mmax):
        new_params = [params[0],params[1],params[2],params[3],m,m]
        sol_pert = pert_sol(new_params)
        mat_1 = K1(sol_pert,sol_backg,new_params)
        #print("K1 = ", mat_1/np.sqrt(2))
        mat += mat_1/np.sqrt(2)
    return mat
