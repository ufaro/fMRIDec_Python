import sys
import numpy as np
import pywt
from math import sqrt
from scipy import linalg
import scipy.io
from nipy.modalities.fmri import hrf, utils
from scipy.linalg import toeplitz, norm



def mad(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)

def soft_thresh(x, l):
    return np.maximum(0,(1-l / (np.maximum(np.abs(x),10 ** -10)))) * x;

def fast_tvprox(z, gamma):

    N = len(z)
    x = np.zeros(N)
    # A:
    k = k0 = kp = km = 0
    vmin = z[0]-gamma
    vmax = z[0]+gamma
    umin = gamma
    umax = -gamma


    while 1:
        # B:
        if(k == N):
            return np.array([vmin+umin])

        # Break condition to avoid overflow...
        if k+1 >= N:
            break

        # C:
        if(z[k+1]+umin < vmin-gamma):
            for i in range(k0, km+1):
                x[i] = vmin
            x[k0] = x[km] = vmin
            k = k0 = km = kp = km+1
            vmin = z[k]
            vmax = z[k]+(2*gamma)
            umin = gamma
            umax = -gamma
        # D:
        elif(z[k+1]+umax > vmax+gamma):
            for i in range(k0, kp+1):
                x[i] = vmax
            x[k0] = x[km] = x[kp] = vmax
            k = k0 = km = kp = kp+1
            vmin = z[k]-(2*gamma)
            vmax = z[k]
            umin = gamma
            umax = -gamma
        # E:
        else:
            k = k+1
            umin = umin +z[k] - vmin
            umax = umax + z[k] - vmax
            # F:
            if(umin >= gamma):
                vmin = vmin + ((umin -gamma)/(k-k0+1))
                umin = gamma
                km = k
            if(umax <= -gamma):
                vmax = vmax+((umax + gamma)/(k-k0+1))
                umax = -gamma
                kp = k
        # G:
        if k >= N:
        # H:
            if(umin < 0):
                for i in range(k0, km+1):
                    x[i] = vmin
                k = k0 = km = km + 1
                vmin = z[k]
                umin = gamma
                umax = z[k] + gamma - vmax
                continue
            # I:
            elif(umax > 0):
                for i in range(k0, kp+1):
                    x[i] = vmax
                k = k0 = kp = kp+1
                vmax = z[k]
                umax = -gamma
                umin = z[k]-gamma-vmin
                continue
            else:
                for i in range(k0, N):
                    x[i] = vmin+(umin/(k-k0+1))
                break

    return x




def GPFM(x,H,Ht,maxeig,lam,Nit,Nh,condition):
    y= np.concatenate((x.T,np.zeros([Nh-1])),axis=0)
    k, t = 1, 1
    s = np.zeros([len(x)])
    u = np.zeros([len(x)])
    if(condition == "blocks"):
      for _ in xrange(Nit):
        u_l = u
        z = Ht.dot(y) / maxeig + s - Ht.dot(H.dot(s)) / maxeig
        u = fast_tvprox(z , lam / maxeig)
        t_l = t
        t = (1. + sqrt(1. +4. * t ** 2)) / 2.
        s = u + ((t_l - 1) / t) * (u - u_l) 
    # E:
    elif(condition == "spikes"):
        tau = lam / maxeig
        for _ in xrange(Nit):
          u_l = u
          z = Ht.dot(y) / maxeig + s - Ht.dot(H.dot(s)) / maxeig
          u = soft_thresh(z , lam / maxeig)
          t_l = t
          t = (1. + sqrt(1. +4. * t ** 2)) / 2.
          s = u + ((t_l - 1) / t) * (u - u_l)

    return u


def volumap(TCN, H, Nh, lam0=1, Nit=100, condition = "blocks"):
    Ht = H.T
    maxeig = norm(H) ** 2  # Lipschitz constant
    nor    = 1 / .6745
    if (len(TCN.shape) == 1):
            y = TCN
            cA, cD = pywt.dwt(y, 'db6')
            lam = mad(cD,0) * nor * lam0
            TC = GPFM(y,H,Ht,maxeig,lam,Nit,Nh,condition)
    elif (len(TCN.shape) == 2):
        T  = TCN.shape[0]
        NbrVoxels = TCN.shape[1] 
        TC = np.zeros_like(TCN)
        for v in range(NbrVoxels):
            y = np.zeros([T])
            y = TCN[:,v]
            cA, cD = pywt.dwt(y, 'db3')
            lam = mad(cD,0) * nor * lam0
            u = GPFM(y,H,Ht,maxeig,lam,Nit,Nh,condition)
            TC[:,v]=u
            if (v % 1000 == 0):
              print(v)

    return TC



def MakeHrfToeplitz(T,dt,du=32):
    t_vec = np.arange(0, du, dt)
    h = hrf.spmt(t_vec)

    hrfunc = np.concatenate((h.T, np.zeros([T-1])),axis=0)
    c      = np.concatenate((np.array([hrfunc[1]]), np.zeros([T-1])),axis=0)
    xConv  = toeplitz(c, hrfunc)
    H      =   xConv.T  
    Nh  = len(h)

    return H, Nh






