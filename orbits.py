"""
Created on Wed Feb 14 21:03:39 2018

@author: Dani Selva
"""

import numpy as np

## Constants
G = 6.674e-11

M_sun = 1.989e30
M_earth = 5.972e24
M_mars = 6.39e23
M_saturn = 5.6834e26
M_titan = 1345.5e20
M_moon = 7.34767309e22
M_venus = 4.8675e24

mu_sun = G * M_sun
#mu_earth = G*M_earth
mu_earth = 3.986004415e14
mu_mars = G * M_mars
mu_saturn = G * M_saturn
mu_titan = G * M_titan
mu_moon = G * M_moon
mu_venus = G * M_venus

R_earth = 6378.1363
R_sun = 695508
R_mars = 3389
R_saturn = 58232
R_titan = 2575
R_moon = 1737
R_venus = 6051.8

AU = 1.496e11
a_earth = 1 * AU
a_mars = 1.524 * AU
a_saturn = 9.557 * AU
a_moon = 0.3844e9
a_venus = 0.723 * AU

#J2_earth = 1.08262668e-3
J2_earth = 1.0826362e-3
J2_venus = 4.458e-6


## Period and utilities

def hkm2a(hkm, R=R_earth):
    return 1000 * (R + hkm)


def a2hkm(a_m, R=R_earth):
    return a_m / 1000 - R


def T2a(T, mu=mu_earth):
    n = 2 * np.pi / T
    a = (mu / (n ** 2)) ** (1.0 / 3)
    return a


def T2h(T, mu=mu_earth, R=R_earth):
    return (T2a(T, mu) - 1000 * R) / 1000


def a2T(a, mu=mu_earth):
    T = 2 * np.pi * np.sqrt(a ** 3 / mu)
    return T


def a2n(a, mu=mu_earth):
    return 2 * np.pi / a2T(a, mu)


def h2T(h_km, mu=mu_earth, R=R_earth):
    a = hkm2a(h_km, R)
    n = np.sqrt(mu / (a ** 3))
    T = 2 * np.pi / n
    return T


def n2a(n, mu=mu_earth):
    a = (mu / n ** 2) ** (1 / 3)
    return a


def n2h(n, mu=mu_earth, R=R_earth):
    a = n2a(n, mu)
    return a / 1000 - R


def h2n(h, mu=mu_earth, R=R_earth):
    return T2n(h2T(h, mu, R))


def n2T(n):
    return 2 * np.pi / n


def T2n(T):
    return 2 * np.pi / T


### Orbital velocity
def orb_vel(a, r, mu=mu_earth):
    v = np.sqrt(mu * (2.0 / r - 1.0 / a))
    return v


def orb_vel_c(h_km, mu=mu_earth, R=R_earth):
    v = np.sqrt(mu / hkm2a(h_km, R))
    return v


def hyperb_vel(r, Vinf, mu=mu_earth):
    V = np.sqrt(Vinf ** 2 + 2 * mu / r)
    return V


def escape_vel(r, mu=mu_earth):
    return np.sqrt(2 * mu / r)


### Hyperbolic orbit utilities
def Vinf2a(Vinf, mu=mu_earth):
    a = -mu / (2 * Vinf ** 2)
    return a


def hyp_e(rp, Vinf, mu=mu_earth):
    return 1 + rp / mu * (Vinf ** 2)


def theta_asymp(e):
    return np.degrees(np.arccos(-1 / e))


## repeat ground track
def Td2a(Td, inc, mu=mu_earth, R=R_earth, J2=J2_earth):
    from perturbations import deltan, SSOh2i, omegadot
    e = 0
    T0 = Td
    n0 = 2 * np.pi / T0
    a0 = T2a(T0, mu)
    h = a0 - R
    da = 1e5
    nn = 1
    #print(a0,da,inc)
    while abs(da) > 1 and nn < 100:
        dn = deltan(a0, inc, e)
        odot = omegadot(a0, inc, e, mu, R, J2)
        dT = -T0 * (odot + dn) / n0
        T0 = Td - dT
        n0 = 2 * np.pi / T0
        a = T2a(T0, mu)
        da = a - a0
        a0 = a
        h = a2hkm(a0, R)
        inc = SSOh2i(h, e, mu, R)
        print(a0, da, inc)
        nn = nn + 1
    return a0


def triplet2h(i, j, k, mu=mu_earth, R=R_earth, J2=J2_earth):
    from perturbations import SSOh2i
    e = 0
    kk = i + j / k
    Td = 86400 / kk
    h0 = T2h(Td, mu, R)
    inc = SSOh2i(h0, e, mu, R)
    a = Td2a(Td, inc, mu, R, J2)
    h = a / 1000 - R
    #print(h)
    return h


def repeat_cycle2h(cycle, h0):
    orb_day = 86400 / h2T(h0)
    i = round(orb_day)
    k = cycle
    j = round(k * (orb_day - i))
    print(i, j, k)
    h = triplet2h(i, j, k)
    return h


def calc_SSO_repeat_orbits():
    from fractions import Fraction
    vec_h = []
    vec_trip = []
    vec_k = []
    vec_j = []
    n = 1
    ok_triplets = {}
    for i in range(14, 17):
        for j in range(-13, 13):
            for k in range(1, 30):
                if j == 0 and k > 1:
                    continue
                if k > abs(j):
                    f = Fraction(j, k)
                    tr = [i, f.numerator, f.denominator]
                    print(tr)
                    if tr not in ok_triplets.values():
                        ok_triplets[n] = tr
                        n = n + 1

    #    return ok_triplets
    for f in ok_triplets.values():
        i = f[0]
        j = f[1]
        k = f[2]
        h = triplet2h(i, j, k)
        if h < 640 and h > 320:
            print(i, j, k, h)
            vec_trip.append([i, j, k])
            vec_k.append(k)
            vec_h.append(h)
            vec_j.append(j)

    return vec_k, vec_trip, vec_h, vec_j


def plot_SSO_repeat_orbits(ks, trs, hs, vec_j):
    import pylab as py
    #plt.scatter(vec_k,vec_h)
    #plt.savefig('repeat.png')
    js = list(map(lambda x: str(x), vec_j))
    for i, s in enumerate(js):
        py.scatter(ks[i], hs[i], marker=r"$ {} $".format(s), s=100, c=[0, 0, 0])
    py.savefig('repeat3.png')
