import Solutions as sol
from numpy import pi

def avg_r(params):
    r_0 , q , mmax = params[0] , params[1] , params[2]
    r_1 = 0
    for m in range(1,mmax):
        r_1 += abs(sol.r_m([r_0,q,m]))/m
    return r_1/pi

def avg_eccentricity(params):
    r_0 = params[0]
    return avg_r(params)*pi/r_0
