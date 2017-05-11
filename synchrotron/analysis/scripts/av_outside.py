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
    print("\nOptimising coefficients below bucket")
    below = np.where(H < 0)
    hb = -H[below]
    cb = coef[below]
    
    cb /= coef.max()
    hb /= H.max()
    
    s = cb.size
    m = cb > 1e-9
    hb = hb[m]
    cb = cb[m]
    print("Data {} -> {} data".format(s, cb.size))
    
    p_init = [-1.0, 1.0, 0]
    res = least_squares(err_f, p_init, loss='linear', f_scale=1e-2, args=(hb, cb), jac='2-point')
    print("----")
    print(res)
    print("----")
    print("Fit: ", *["{:.3f},".format(v) for v in res.x])
    
    fit = lambda x : model(x, *res.x)
    plot_coefficients(hb, cb, plot_type="scale", curve=fit, info="below")
    return (hb, cb)

def above():
    print("\nOptimising coefficients above bucket")
    above = np.where(H >= 0)
    ha = H[above]
    ca = coef[above]
    
    ca /= coef.max()
    ha /= H.max()
    
    s = ca.size
    m = ca > 1e-5
    m *= ha > 1e-6
    m *= np.invert((ca > 0.1)*(ha > 0.25))
    # m *= np.invert((ca < 3e-2)*(ha < 0.5))
    # m *= np.invert((ha > 0.8)*(ca > 0.1))
    ha = ha[m]
    ca = ca[m]
    print("Data: {} -> {}".format(s, ca.size))
    
    p_init = [-1.0, 1.0, 0]
    res = least_squares(err_f, p_init, max_nfev=1000, loss='cauchy', f_scale=1e-2, args=(ha, ca), jac='2-point')
    print("----")
    print(res)
    print("----")
    print("Fit: ", *["{:.3f},".format(v) for v in res.x])
    
    fit = lambda x : model(x, *res.x)
    plot_coefficients(ha, ca, plot_type="scale", curve=fit, info="above")
    return (ha, ca)

ha, ca = above()
hb, cb = below()

h = list(hb) + list(ha)
c = list(cb) + list(ca)

r = sum(ca)/(sum(ca) + sum(cb))
print("Ratio particles\n\tbelow: {:.4f}\n\tabove: {:.4f}".format(1 - r, r))

plot_coefficients(H, coef, plot_type='bar', info="Initial", block=False)
plot_coefficients(h, c, plot_type='bar', info="Pruned")
