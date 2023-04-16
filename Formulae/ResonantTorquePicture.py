import Potentials as pt
import numpy as np
from scipy.optimize import bisect

def dDdlnr(m,n):
    return -3*n*n/(m+1)

def LR_location(m,n):
    return ((m+1)/n)**(2./3.)

def Psi(params):
    r, m, n = params[0], params[4], params[5]
    phi = pt.compute_eccentric_potential(params)
    dphi = pt.compute_eccentric_potential_d1(params)
    omega = r**(-1.5)
    omega_p = n/m
    return r*dphi + 2*omega*phi/(omega - omega_p)
    
def Torque(params):
    m, n = params[4], params[5]
    psi = Psi(params)
    dD_dlnr = dDdlnr(m,n)
    torque = -m*np.pi*np.pi*psi*psi/dD_dlnr
    return torque

def zeta_deltat(fluid_params):
    params, alpha, chi = fluid_params[0], fluid_params[1], fluid_params[2]
    Tbar_mn = Torque(params)
    r = params[0]
    return Tbar_mn/r/alpha/chi/chi

def zeta_T(fluid_params):
    zeta_dt = zeta_deltat(fluid_params)
    return zeta_dt/3./np.pi

def unpert_sol(params):
    r0 = params[0]
    return np.array([r0,r0**(-1.5),0,np.sqrt(r0)])

def pert_sol(params):
    m , n = params[4] , params[5]
    phi = pt.compute_eccentric_potential(params)
    dphi = pt.compute_eccentric_potential_d1(params)
    backg_sol = unpert_sol(params)
    r0 , omega0 = backg_sol[0] , backg_sol[1]
    D = (omega0)**2 - (m*omega0 - n)**2
    return [-(2.*phi*(1 - n/(m*omega0))/r0 + dphi)/np.abs(D),-m*phi/(m*omega0 - n)]

def root_func(params):
    m,n = params[4],params[5]
    rLR = LR_location(m,n)
    amp = pert_sol(params)[0]
    return np.log10( abs(amp/(params[0] - rLR)) )

def find_rgap(params,guess_vars):
    start,fin = guess_vars[0],guess_vars[1]
    rfunc = lambda x: root_func([x,params[1],params[2],params[3],params[4],params[5]])
    rgap = bisect(rfunc,start,fin)
    return rgap
