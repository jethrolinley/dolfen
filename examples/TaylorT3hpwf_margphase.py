
from astropy import constants as const
import numpy as np


# initialise a few constants
pi=np.pi
pi2=pi**2
SWSHpf = (1./8.)*np.sqrt(5./pi)  # some prefactors ('pf') to prevent recalculation
H22pf = 8.*np.sqrt(np.pi/5.) 
G = const.G.value
c = const.c.value
M_sun = const.M_sun.value
mpc = const.pc.value*1.0e6
G_c3 = G/c**3
mpc_c = mpc/c

def TaylorT3hpwaveform_margphase(t,geocent_time,mass_1,mass_2,a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl, theta_jn, luminosity_distance, psi, ra, dec, phase=0.0):
    # iota = theta_jn when aligned spins!

    # h plus waveform only...
    # 10 params for aligned spin BHBs
    # set scaled params eta, M, and so on
    # time series, coalescence time, mass1,mass2,spin2,spin2,distance,coalescence phase,inclination,azimuth,sky position(x1,x2),sky orientatin(x3)
    # from https://arxiv.org/pdf/2004.08302.pdf
    # https://arxiv.org/pdf/0802.1249.pdf
    # https://arxiv.org/pdf/gr-qc/0105099.pdf
    # https://journalaporg/prd/abstract/10.1103/PhysRevD.80.084043
    # https://journalaporg/prd/abstract/10.1103/PhysRevD.87.044009

    # say units of params in comments!!

    M = mass_1+mass_2
    Mtot = M*M_sun*G_c3
    mass_1 /= M
    mass_2 /= M
    M = 1.
    luminosity_distance *= mpc_c
    cosi = np.cos(theta_jn)
    n = mass_1*mass_2
    dm = np.sqrt(1.-4.*n)   #'delta m'

    tc_s_t = (geocent_time-t)/Mtot
    tau = tc_s_t*n/5
    tau1_8 = tau**0.125
    tau2_8 = tau1_8*tau1_8
    tau3_8 = tau1_8*tau2_8
    tau5_8 = tau3_8*tau2_8
    theta = 1./tau1_8
    theta2 = theta*theta
    theta3 = theta*theta2
    #theta4 = theta2*theta2
    theta5 = theta2*theta3
    theta6 = theta3*theta3
    theta7 = theta5*theta2
    theta8 = theta5*theta3
    theta9 = theta6*theta3
    theta10 = theta5*theta5

    xsxx       = a_1-a_2
    xpxx       = a_1+a_2
    mx         = mass_1*a_1
    mmxx       = mass_2*a_2
    mxpmmxx    = mx+mmxx
    nxpxx      = n*xpxx
    x2         = a_1*a_1
    xx2        = a_2*a_2
    xxx        = a_1*a_2
    mx2pmmxx2  = mass_1*x2 + mass_2*xx2
    n2         = n*n
    npi        = n*pi
    ndxsxx     = n*dm*xsxx
    n2xpxx     = n*n*xpxx
    y          = np.euler_gamma
    npi2       = n*pi2
    n2pi       = n2*pi
    n3         = n2*n

    nx2pxx2    = n*(x2+xx2)
    nxxx       = n*xxx
    ndx2sxx2   = n*dm*(x2-xx2)
    n2x2pxx2   = n*nx2pxx2
    n2xxx      = n*nxxx
    mx3pmmxx3  = mx*x2 + mmxx*xx2
    npix2pxx2  = npi*(x2+xx2)
    npixxx     = npi*xxx
    nx3pxx3    = n*(x2*a_1 + xx2*a_2)
    ndx3sxx3   = n*dm*(x2*a_1 - xx2*a_2)
    n2x3pxx3   = n*nx3pxx3
    dxsxx      = dm*xsxx
    n2dxsxx    = n2*dxsxx

    #wh0 = 1.
    #wh1 = 0.
    wh2 = 743/2688 + 11*n/32
    wh3 = -3*pi/10 + 113/160*mxpmmxx - (19/80)*nxpxx
    wh4 = 1855099/14450688 - 243*mx2pmmxx2/1024 + 56975*n/258048 + (3/1024)*n*(81*x2-158*xxx+81*xx2) + 371*n2/2048
    wh5 = -7729*pi/21504 + 146597*mxpmmxx/64512 + 13*npi/256 - 1213*nxpxx/1152 + 7*ndxsxx/128 - 17*n2xpxx/128
    wh6 = -720817631400877/288412611379200 + 107*y/280 + 53*pi2/200 - 6127*pi*mxpmmxx/6400 - 16928263*mx2pmmxx2/68812800 \
            + 25302017977*n/4161798144 - 451*npi2/2048 + 1051*pi*nxpxx/3200 + 23281001*nx2pxx2/68812800 - 377345*nxxx/1376256 + 453767*ndx2sxx2/4915200 \
            - 30913*n2/1835008 + 335129*n2x2pxx2/2457600 - 488071*n2xxx/1228800 + (107/280)*np.log(2*theta)
    wh7 = -188516689*pi/433520640 + 6579635551*mxpmmxx/650280960 + 3663*pi*mx2pmmxx2/5120 - 67493*mx3pmmxx3/81920 \
            -97765*npi/258048 - 1496368361**nxpxx/185794560 - 3663*npix2pxx2/5120  + 3537*npixxx/2560 + 206917*nx3pxx3/163840 \
            -192709*n*xxx*mxpmmxx/81920 - 28633921*ndxsxx/12386304 + 71931*ndx3sxx3/163840 + 141769*n2pi/1290240 \
            -840149*n2xpxx/15482880 - 2219*n2x3pxx3/40960 + 1343*n2xpxx*xxx/40960 + 2591*n2dxsxx/46080 - 12029*n*n2xpxx/92160


    omega22 = theta3 + wh2*theta5 + wh3*theta6 + wh4*theta7 + wh5*theta8 + wh6*theta9 + wh7*theta10
    psi22 = -tau5_8  - 5.*(wh2*tau3_8/3. + wh4*tau1_8 - wh6*theta) - 2.5*(wh3*tau2_8 - wh7*theta2) - 0.625*wh5*np.log(tc_s_t)
    psi22 /= n
    psi22 += phase

    x1_2 = (0.125*omega22)**(1./3.)
    x    = x1_2*x1_2

    hh2 = -107/42 + 55*n/42
    hh3 = 2*pi - 2*xpxx/3 + 2*dxsxx/(3*M) + 2*nxpxx/3
    hh4 = -2173/1512 - 1069*n/216 + 2047*n2/1512 + mx2pmmxx2 - n*xsxx**2

    H22 = 1. + x*(hh2 + x1_2*hh3 + x*hh4)
    H22 *= H22*n*x
    H22 *= Mtot/luminosity_distance
    h22cos = H22*np.exp(2.j*psi22)
    Y22 = (1.-cosi)**2*np.exp(-2.j*psi)
    Y2_2 = (1.+cosi)**2*np.exp(2.j*psi)
    h_cos = np.real(h22cos*Y22 + np.conj(h22cos)*Y2_2)

    h22sin = H22*np.exp(2.j*(psi22-np.pi/4.))
    h_sin = np.real(h22sin*Y22 + np.conj(h22sin)*Y2_2)

    return h_cos, h_sin



def TaylorT3hpwaveform(t,geocent_time,mass_1,mass_2,a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl, theta_jn, luminosity_distance, phase, psi, ra, dec):#, zeta, EdGB):

    M = mass_1+mass_2
    Mtot = M*M_sun*G_c3
    mass_1 /= M
    mass_2 /= M
    M = 1.
    luminosity_distance *= mpc_c
    cosi = np.cos(theta_jn)
    n = mass_1*mass_2#/M**2          #'eta'
    dm = np.sqrt(1.-4.*n)   #'delta m'

    tc_s_t = (geocent_time-t)/Mtot
    tau = tc_s_t*n/5
    tau1_8 = tau**0.125
    tau2_8 = tau1_8*tau1_8
    tau3_8 = tau1_8*tau2_8
    tau5_8 = tau3_8*tau2_8
    theta = 1./tau1_8
    theta2 = theta*theta
    theta3 = theta*theta2
    #theta4 = theta2*theta2
    theta5 = theta2*theta3
    theta6 = theta3*theta3
    theta7 = theta5*theta2
    theta8 = theta5*theta3
    theta9 = theta6*theta3
    theta10 = theta5*theta5

    xsxx       = a_1-a_2
    xpxx       = a_1+a_2
    mx         = mass_1*a_1
    mmxx       = mass_2*a_2
    mxpmmxx    = mx+mmxx
    nxpxx      = n*xpxx
    x2         = a_1*a_1
    xx2        = a_2*a_2
    xxx        = a_1*a_2
    mx2pmmxx2  = mass_1*x2 + mass_2*xx2
    n2         = n*n
    npi        = n*pi
    ndxsxx     = n*dm*xsxx
    n2xpxx     = n*n*xpxx
    y          = np.euler_gamma
    npi2       = n*pi2
    n2pi       = n2*pi
    n3         = n2*n

    nx2pxx2    = n*(x2+xx2)
    nxxx       = n*xxx
    ndx2sxx2   = n*dm*(x2-xx2)
    n2x2pxx2   = n*nx2pxx2
    n2xxx      = n*nxxx
    mx3pmmxx3  = mx*x2 + mmxx*xx2
    npix2pxx2  = npi*(x2+xx2)
    npixxx     = npi*xxx
    nx3pxx3    = n*(x2*a_1 + xx2*a_2)
    ndx3sxx3   = n*dm*(x2*a_1 - xx2*a_2)
    n2x3pxx3   = n*nx3pxx3
    dxsxx      = dm*xsxx
    n2dxsxx    = n2*dxsxx

    #wh0 = 1.
    #wh1 = 0.
    wh2 = 743/2688 + 11*n/32
    wh3 = -3*pi/10 + 113/160*mxpmmxx - (19/80)*nxpxx
    wh4 = 1855099/14450688 - 243*mx2pmmxx2/1024 + 56975*n/258048 + (3/1024)*n*(81*x2-158*xxx+81*xx2) + 371*n2/2048
    wh5 = -7729*pi/21504 + 146597*mxpmmxx/64512 + 13*npi/256 - 1213*nxpxx/1152 + 7*ndxsxx/128 - 17*n2xpxx/128
    wh6 = -720817631400877/288412611379200 + 107*y/280 + 53*pi2/200 - 6127*pi*mxpmmxx/6400 - 16928263*mx2pmmxx2/68812800 \
            + 25302017977*n/4161798144 - 451*npi2/2048 + 1051*pi*nxpxx/3200 + 23281001*nx2pxx2/68812800 - 377345*nxxx/1376256 + 453767*ndx2sxx2/4915200 \
            - 30913*n2/1835008 + 335129*n2x2pxx2/2457600 - 488071*n2xxx/1228800 + (107/280)*np.log(2*theta)
    wh7 = -188516689*pi/433520640 + 6579635551*mxpmmxx/650280960 + 3663*pi*mx2pmmxx2/5120 - 67493*mx3pmmxx3/81920 \
            -97765*npi/258048 - 1496368361**nxpxx/185794560 - 3663*npix2pxx2/5120  + 3537*npixxx/2560 + 206917*nx3pxx3/163840 \
            -192709*n*xxx*mxpmmxx/81920 - 28633921*ndxsxx/12386304 + 71931*ndx3sxx3/163840 + 141769*n2pi/1290240 \
            -840149*n2xpxx/15482880 - 2219*n2x3pxx3/40960 + 1343*n2xpxx*xxx/40960 + 2591*n2dxsxx/46080 - 12029*n*n2xpxx/92160



    omega22 = theta3 + wh2*theta5 + wh3*theta6 + wh4*theta7 + wh5*theta8 + wh6*theta9 + wh7*theta10
    psi22 = -tau5_8  - 5.*(wh2*tau3_8/3. + wh4*tau1_8 - wh6*theta) - 2.5*(wh3*tau2_8 - wh7*theta2) - 0.625*wh5*np.log(tc_s_t)
    psi22 /= n
    psi22 += phase

    omega22 = theta3 + wh2*theta5 + wh3*theta6 + wh4*theta7 + wh5*theta8 + wh6*theta9 + wh7*theta10
    x1_2 = (0.125*omega22)**(1./3.)
    x    = x1_2*x1_2

    hh2 = -107/42 + 55*n/42
    hh3 = 2*pi - 2*xpxx/3 + 2*dxsxx/(3*M) + 2*nxpxx/3
    hh4 = -2173/1512 - 1069*n/216 + 2047*n2/1512 + mx2pmmxx2 - n*xsxx**2

    H22 = 1. + x*(hh2 + x1_2*hh3 + x*hh4)
    H22 *= H22*n*x
    h22 = H22*np.exp(2.j*psi22)
    Y22 = (1.-cosi)**2*np.exp(-2.j*psi)
    Y2_2 = (1.+cosi)**2*np.exp(2.j*psi)

    h = h22*Y22 + np.conj(h22)*Y2_2
    h *= Mtot/luminosity_distance

    return np.real(h) 


def TaylorT3hpwaveformMaxFreq(t,geocent_time,mass_1,mass_2,a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl, theta_jn, luminosity_distance, phase, psi, ra, dec, **kwargs):
    M = mass_1+mass_2
    Mtot = M*M_sun*G_c3
    mass_1 /= M
    mass_2 /= M
    M = 1.
    luminosity_distance *= mpc_c
    cosi = np.cos(theta_jn)
    n = mass_1*mass_2/M**2          #'eta'
    dm = np.sqrt(1.-4.*n)   #'delta m'

    tc_s_t = (geocent_time-t)/Mtot
    tau = tc_s_t*n/5
    tau1_8 = tau**0.125
    tau2_8 = tau1_8*tau1_8
    tau3_8 = tau1_8*tau2_8
    tau5_8 = tau3_8*tau2_8
    theta = 1./tau1_8
    theta2 = theta*theta
    theta3 = theta*theta2
    #theta4 = theta2*theta2
    theta5 = theta2*theta3
    theta6 = theta3*theta3
    theta7 = theta5*theta2
    theta8 = theta5*theta3
    theta9 = theta6*theta3
    theta10 = theta5*theta5

    xsxx       = a_1-a_2
    xpxx       = a_1+a_2
    mx         = mass_1*a_1
    mmxx       = mass_2*a_2
    mxpmmxx    = mx+mmxx
    nxpxx      = n*xpxx
    x2         = a_1*a_1
    xx2        = a_2*a_2
    xxx        = a_1*a_2
    mx2pmmxx2  = mass_1*x2 + mass_2*xx2
    n2         = n*n
    npi        = n*pi
    ndxsxx     = n*dm*xsxx
    n2xpxx     = n*n*xpxx
    y          = np.euler_gamma
    npi2       = n*pi2
    n2pi       = n2*pi
    n3         = n2*n

    nx2pxx2    = n*(x2+xx2)
    nxxx       = n*xxx
    ndx2sxx2   = n*dm*(x2-xx2)
    n2x2pxx2   = n*nx2pxx2
    n2xxx      = n*nxxx
    mx3pmmxx3  = mx*x2 + mmxx*xx2
    npix2pxx2  = npi*(x2+xx2)
    npixxx     = npi*xxx
    nx3pxx3    = n*(x2*a_1 + xx2*a_2)
    ndx3sxx3   = n*dm*(x2*a_1 - xx2*a_2)
    n2x3pxx3   = n*nx3pxx3
    dxsxx      = dm*xsxx
    n2dxsxx    = n2*dxsxx

    #wh0 = 1.
    #wh1 = 0.
    wh2 = 743/2688 + 11*n/32
    wh3 = -3*pi/10 + 113/160*mxpmmxx - (19/80)*nxpxx
    wh4 = 1855099/14450688 - 243*mx2pmmxx2/1024 + 56975*n/258048 + (3/1024)*n*(81*x2-158*xxx+81*xx2) + 371*n2/2048
    wh5 = -7729*pi/21504 + 146597*mxpmmxx/64512 + 13*npi/256 - 1213*nxpxx/1152 + 7*ndxsxx/128 - 17*n2xpxx/128
    wh6 = -720817631400877/288412611379200 + 107*y/280 + 53*pi2/200 - 6127*pi*mxpmmxx/6400 - 16928263*mx2pmmxx2/68812800 \
            + 25302017977*n/4161798144 - 451*npi2/2048 + 1051*pi*nxpxx/3200 + 23281001*nx2pxx2/68812800 - 377345*nxxx/1376256 + 453767*ndx2sxx2/4915200 \
            - 30913*n2/1835008 + 335129*n2x2pxx2/2457600 - 488071*n2xxx/1228800 + (107/280)*np.log(2*theta)
    wh7 = -188516689*pi/433520640 + 6579635551*mxpmmxx/650280960 + 3663*pi*mx2pmmxx2/5120 - 67493*mx3pmmxx3/81920 \
            -97765*npi/258048 - 1496368361**nxpxx/185794560 - 3663*npix2pxx2/5120  + 3537*npixxx/2560 + 206917*nx3pxx3/163840 \
            -192709*n*xxx*mxpmmxx/81920 - 28633921*ndxsxx/12386304 + 71931*ndx3sxx3/163840 + 141769*n2pi/1290240 \
            -840149*n2xpxx/15482880 - 2219*n2x3pxx3/40960 + 1343*n2xpxx*xxx/40960 + 2591*n2dxsxx/46080 - 12029*n*n2xpxx/92160

    # phase, (2,2) mode
    # orbital phase!!
    psi22 = -tau5_8  - 5.*(wh2*tau3_8/3. + wh4*tau1_8 - wh6*theta) - 2.5*(wh3*tau2_8 - wh7*theta2) - 0.625*wh5*np.log(tc_s_t)
    psi22 /= n
    psi22 += phase

    omega22 = theta3 + wh2*theta5 + wh3*theta6 + wh4*theta7 + wh5*theta8 + wh6*theta9 + wh7*theta10
    return np.max(np.abs(omega22))/(8*np.pi*Mtot)



def T3NyquistFrequency(times,**kwargs):
    return 2.*TaylorT3hpwaveformMaxFreq(times,**kwargs)


