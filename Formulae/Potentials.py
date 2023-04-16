import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial,lpmv

## params get input into this module as
## r, q, i, ecc, m, N
## params get parsed into potentials as
## r, q, i, ecc, m, mu, n


## Define the Wigner little-d matrix

def wignerd(j,mp,m,beta):
  prefactor = np.sqrt(factorial(j+mp)*factorial(j-mp)*factorial(j+m)*factorial(j-m))
  temp = 0
  smin = max(0,m-mp)
  smax = min(j+m,j-mp)
  th = 0.5*beta
  cth = np.cos(th)
  sth = np.sin(th)
  s = smin
  while s < smax+1:
    f1 , f2 , f3 = j + m - s , mp - m + s , j - mp - s
    if f1 < 0 or f2 < 0 or f3 < 0:
      break
    else:
      temp1 = (-1)**(mp - m + s)
      temp2 = cth**(2*j + m - mp - 2*s)
      temp3 = sth**(mp - m + 2*s)
      temp += temp1*temp2*temp3/factorial(f1)/factorial(f2)/factorial(f3)/factorial(s)
      s += 1
  return prefactor * temp

## define the W coefficient from integrating the Legendre Polynomial

def W(l,m):
    return W1(l,m)

def W1(l,m):
  return ( ( factorial(l-m)/factorial(l+m) ) **0.5) * lpmv(m,l,0)

def W2(l,m):
    lpm = (l+m)
    lmm = (l-m)
    lpmh = lpm*0.5
    lmmh = lmm*0.5
    fp = factorial(lpm)
    fph = factorial(lpmh)
    fm = factorial(lmm)
    fmh = factorial(lmmh)
    if lpm%2 == 0:
        return ((-1)**lpmh)*np.sqrt(fp*fm)/((2**l)*fph*fmh)
    return 0

#for m in range(2,10):
#    for l in range(m,m+10):
#        print("(l,m):",l,",",m)
#        print("W with the legendre polynomial:", W1(l,m))
#        print("W with the direct method:", W2(l,m))
        
        

## define the eccentricity corrections

def Clmn(l,m,n,ecc):
    ecc2 = ecc*ecc
    ecc3 = ecc2*ecc
    ecc4 = ecc3*ecc
    l2 = l*l
    l3 = l2*l
    l4 = l3*l
    m2 = m*m
    m3 = m2*m
    if n == 0:
        return 1
    if n == 1:
        return (1/2)*ecc*(-l + 2*m)
    if n == -1:
        return (1/2)*ecc*(-l-2*m)
    if n == 2:
        return (1/8)*ecc2*(l2 - (3 + 4*m)*l + (4*m + 5)*m)
    if n == -2:
        return (1/8)*ecc2*(l2 - (3 - 4*m)*l + (4*m - 5)*m)
    if n == 3:
        return (1/48)*ecc3*(-l3 + 3*(3 + 2*m)*l2 - (12*m2 + 33*m + 17)*l + 2*m*(4*m2 + 15*m + 13))
    if n == -3:
        return (1/48)*ecc3*(-l3 + 3*(3 - 2*m)*l2 - (12*m2 - 33*m + 17)*l + 2*m*(-4*m2 + 15*m - 13))
    if n == 4:
        return (1/384)*ecc4*(l4 - 2*(9 + 4*m)*l3 + (24*m2 + 102*m + 95)*l2 - 2*l*(16*m3 + 96*m2 + 165*m + 71) + (16*m3 + 120*m2 + 283*m + 206)*m)
    if n == -4:
        return (1/384)*ecc4*(l4 - 2*(9 - 4*m)*l3 + (24*m2 - 102*m + 95)*l2 - 2*l*(-16*m3 + 96*m2 - 165*m + 71) + (16*m3 - 120*m2 + 283*m - 206)*m)
    return "bad choice of params"

## define the selection rules for a choice of params
## Although there are an infinite choice of values, we will use only use
## values for mu , n such that n is in the range [-4,4]

def select_rules(params):
    N = params[5]
    triplet_params = []
    for n in range(-4,5):
        mu = N - n
        new_params = [params[0],params[1],params[2],params[3],params[4],mu,n]
        triplet_params.append(new_params)
    return triplet_params

def compute_eccentric_potential(params):
    triplet_params = select_rules(params)
    phi_mN = 0
    for new_params in triplet_params:
        phi_mN += modewise_Phi_grav(new_params)
        #print("%d,%d,%d,%.3e"%(new_params[4],new_params[5],new_params[6],phi_mN))
    return phi_mN

def compute_eccentric_potential_d1(params):
    triplet_params = select_rules(params)
    dphi_mN = 0
    for new_params in triplet_params:
        dphi_mN += modewise_dPhi_grav(new_params)
    return dphi_mN

def compute_eccentric_potential_d2(params):
    triplet_params = select_rules(params)
    ddphi_mN = 0
    for new_params in triplet_params:
        ddphi_mN += modewise_ddPhi_grav(new_params)
    return ddphi_mN

## Define the background potential

def backg_Phi(params):
    return -1/params[0]

## Define the background l-polar potential

def backg_multipole_Phi(l,params):
    r , q = params[0] , params[1]
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    W_lm = W(l,0)
    return -2*Q_l*W_lm*W_lm/(r**(l+1))


## define the modewise perturbed potential coefficient and its radial derivative

def modewise_Phi_grav(params):
  r, q, i, ecc , m , mu , n = params[0] , params[1] , params[2] , params[3] , params[4], params[5], params[6]
  lmin = max(2,m,abs(mu))
  lmax = lmin + 6
  #lmax = 2
  pmmun = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*((-q)**l) + q*((1-q)**l))
    wigner = wignerd(l,mu,m,i)
    ecc_term = Clmn(l,mu,n,ecc)
    pmmun -= 2*Q_l*W(l,m)*W(l,mu)*wigner*ecc_term/(r**(l+1))
  return pmmun

def modewise_dPhi_grav(params):
  r , q , i , ecc , m , mu , n = params[0] , params[1] , params[2] , params[3] , params[4] , params[5] , params[6]
  lmin = max(2,m,abs(mu))
  lmax = lmin + 6
  #lmax = 2
  dpmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    wigner = wignerd(l,mu,m,i)
    ecc_term = Clmn(l,mu,n,ecc)
    dpmn += 2*(l+1)*Q_l*W(l,m)*W(l,mu)*wigner*ecc_term/(r**(l+2))
  return dpmn

def modewise_ddPhi_grav(params):
  r , q , i , ecc , m , mu , n = params[0] , params[1] , params[2] , params[3] , params[4] , params[5] , params[6]
  lmin = max(2,m,abs(mu))
  lmax = lmin + 6
  #lmax = 2
  ddpmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    wigner = wignerd(l,mu,m,i)
    ecc_term = Clmn(l,mu,n,ecc)
    ddpmn -= 2*(l+1)*(l+2)*Q_l*W(l,m)*W(l,mu)*wigner*ecc_term/(r**(l+3))
  return ddpmn    

## define the perturbed potential

def pert_Phi(phi,params):
    m = params[4]
    phi_m = modewise_Phi_grav(params)
    return phi_m*np.cos(m*phi)