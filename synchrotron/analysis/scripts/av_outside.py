# NOTE: need to call this file from its local directory for it to work with the import
import sys
sys.path.insert(0,'..')
from lhccomp import *

coef, H = np.loadtxt("av_outside.txt", skiprows=1, usecols=(1, 2), unpack=True)
H /= max(H.max(), -H.min())
coef /= coef.max()

def below():
    print("\nOptimising coefficients below bucket")
    below = np.where(H < 0)
    h = -H[below]
    c = coef[below]
    
    s = c.size
    m = c > 1e-9
    m *= np.invert((h > 0.2)*(c > 0.15))
    m *= np.invert((h < 0.3)*(c < 5e-3))
    h = h[m]
    c = c[m]
    print("Data {} -> {} data".format(s, c.size))

    model = lambda x, *p: np.exp(p[0]*np.power(x, p[1]) + p[2])
    err_f = lambda p, x, y: (model(x, *p) - y)**2
    
    p_init = [-1.0, 1.0, 0]
    res = least_squares(err_f, p_init, loss='linear', f_scale=1e-2, args=(h, c), jac='2-point')
    if not res.success:
        print("----")
        print(res)
        print("----")
    print("Below fit: ", *["{:.3f},".format(v) for v in res.x])
    
    fit = lambda x : model(x, *res.x)
    plot_coefficients(h, c, plot_type="scale", curve=fit, info="below", block=True)
    return (h, c)

def above():
    print("\nOptimising coefficients above bucket")
    above = np.where(H >= 0)
    h = H[above]
    c = coef[above]
    
    s = c.size
    # m = c > 1e-5
    m = c > 3e-3
    m *= h > 1e-6
    m *= np.invert((c > 0.1)*(h > 0.25))
    # m *= np.invert((h > 0.8)*(c > 0.1))
    h = h[m]
    c = c[m]
    print("Data: {} -> {}".format(s, c.size))
    
    model = lambda x, *p: np.exp(p[0]*np.power(x, p[1]) + p[2])
    err_f = lambda p, x, y: (model(x, *p) - y)**2

    p_init = [-1.0, 1.0, 0]
    res = least_squares(err_f, p_init, max_nfev=1000, loss='cauchy', f_scale=1e-2, args=(h, c), jac='2-point')
    if not res.success:
        print("----")
        print(res)
        print("----")
    print("Above fit: ", *["{:.3f},".format(v) for v in res.x])
    
    fit = lambda x : model(x, *res.x)
    plot_coefficients(h, c, plot_type="scale", curve=fit, info="above", block=False)
    return (h, c)

ha, ca = above()
hb, cb = below()

h = list(hb) + list(ha)
c = list(cb) + list(ca)

r = sum(ca)/(sum(ca) + sum(cb))
print("Ratio particles\n\tabove: {:.4f}\n\tbelow: {:.4f}".format(r, 1 - r))

plot_coefficients(H, coef, plot_type='bar', info="Initial", block=False)
plot_coefficients(h, c, plot_type='bar', info="Pruned")
