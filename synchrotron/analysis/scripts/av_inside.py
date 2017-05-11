import sys
sys.path.insert(0, "..")
from lhccomp import *

coef, H = np.loadtxt("av_inside.txt", skiprows=1, usecols=(1, 2), unpack=True)

# h = -H/H.min() # negative values
h = H
c = coef/coef.max()

def model(x, *p):
    return p[0]*x + p[1]

def err(p, x, y):
    return (model(x, *p) - y)**2

m = c > 1e-1
m *= h < -3500
h = h[m]
c = c[m]

res = least_squares(err, [0.0, 0.0], loss='linear', args=(h, c), jac='2-point')
if not res.success:
    print(res)
else:
    print("message:", res.message)
    print("optimality: ", res.optimality)
    print("success:", res.success)
    print("Fit: ", *res.x)

fit = lambda x : model(x, *res.x)
plot_coefficients(h, c, plot_type="scatter", curve=fit)
plot_coefficients(h, c, plot_type='bar')
