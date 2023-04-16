import Potentials as pt
import RotatingSolutions as sol
from numpy import pi

## Define the general kinetic term

def kinetic_energy(r,pr,l):
    return (pr*pr + l*l/r/r)/2 - l

## Define the background Jacobi constant

def backg_C_j(params):
    phase_space = sol.backg_sol_rotating(0,params)
    phi_grav = pt.backg_Phi(params)
    r , pr , l = phase_space[0] , phase_space[1] , phase_space[2]
    return -(kinetic_energy(r,pr,l) + phi_grav)

## Define the background quadrupole Jacobi constant

def backg_quad_C_j(params):
    phase_space = sol.backg_sol_rotating(0,params)
    phi_grav = pt.backg_Phi(params) + pt.backg_multipole_Phi(2,params)
    r , pr , l = phase_space[0] , phase_space[1] , phase_space[2]
    return -(kinetic_energy(r,pr,l) + phi_grav)

## Define the modewise O(\epsilon) corrections to the Jacobi constant

def modewise_C_j(params):
    m = params[4]
    backg = sol.backg_sol_rotating(0,params)
    r0 , l0 = backg[0] , backg[2]
    r02 = r0*r0
    r03 = r02*r0
    r1 = sol.r_m_rotating(params)
    l1 = sol.l_m_rotating(params)
    return -(2/m/pi)*(-l1 + l0*l1/(r02) - l0*l0*r1/(r03) + pt.modewise_Phi_grav(params))

## Define the total O(\epsilon) corrections to the Jacobi constant

def pert_C_j(params):
    mmax = params[4]
    c_j = backg_C_j(params)
    for m in range(1,mmax):
        c_j += modewise_C_j([params[0],params[1],params[2],params[3],m,m])
    return c_j
