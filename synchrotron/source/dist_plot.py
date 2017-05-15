import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("d.txt")
fig, ax = plt.subplots()
bins = np.arange(d.min(), d.max(), (d.max() - d.min())/200)
ax.hist(d, edgecolor='white', bins=bins)
plt.show()
