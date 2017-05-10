# NOTE: need to call this file from its local directory for it to work with the import
import sys
sys.path.insert(0,'..')
from lhccomp import *

def model(x, *p):
    return np.exp(p[0]*np.power(x, p[1]) + p[2])

def err_f(p, x, y):
    return (model(x, *p) - y)**2

coef, H = np.loadtxt("av_outside.txt", skiprows=1, usecols=(1, 2), unpack=True)

def below():
    print("Optimising coefficients below bucket\n")
    below = np.where(H < 0)
    hb = -H[below]
    cb = coef[below]
    
    cb /= cb.max()
    hb /= hb.max()
    
    s = cb.size
    m = np.invert(cb < 1e-9)
    m *= np.invert((hb < 0.3)*(cb < 5e-3))
    m *= np.invert((hb < 0.2)*(cb < 1e-1))
    m *= np.invert((hb < 0.3)*(cb < 2e-2))
    m *= np.invert((hb < 0.1)*(cb < 0.25))
    hb = hb[m]
    cb = cb[m]
    print("Data {} -> {} data".format(s, cb.size))
    
    p_init = [-1.0, 1.0, 0]
    res = least_squares(err_f, p_init, loss='cauchy', f_scale=1e-2, args=(hb, cb), jac='2-point')
    print("----")
    print(res)
    print("----")
    print("Fit: ", *res.x)
    
    fit = lambda x : model(x, *res.x)
    plot_coefficients(hb, cb, plot_type="scale", curve=fit, info="below")
    return (hb, cb)

def above():
    print("Optimising coefficients below bucket\n")
    above = np.where(H >= 0)
    ha = H[above]
    ca = coef[above]
    
    ca /= ca.max()
    ha /= ha.max()
    
    s = ca.size
    m = np.invert(ca < 1e-2)
    m *= ha > 1e-6
    m *= np.invert((ca < 3e-2)*(ha < 0.5))
    m *= np.invert((ha > 0.8)*(ca > 0.1))
    ha = ha[m]
    ca = ca[m]
    print("Data: {} -> {}".format(s, ca.size))
    
    p_init = [-1.0, 1.0, 0]
    res = least_squares(err_f, p_init, loss='cauchy', f_scale=1e-2, args=(ha, ca), jac='2-point')
    print("----")
    print(res)
    print("----")
    print("Fit: ", *res.x)
    
    fit = lambda x : model(x, *res.x)
    plot_coefficients(ha, ca, plot_type="scale", curve=fit, info="above")
    return (ha, ca)

hb, cb = below()
ha, ca = above()

h = list(hb) + list(ha)
c = list(cb) + list(ca)
plot_coefficients(H, coef, plot_type='bar', info="Initial", block=False)
plot_coefficients(h, c, plot_type='bar', info="Pruned")
