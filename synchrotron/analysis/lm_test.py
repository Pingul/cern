import matplotlib.pyplot as plt
import lossmap as lm
from phasespace import PhaseSpace
from settings import settings
import numpy as np

hits = lm.hits_from_collfile(settings.COLL_PATH)
old_lm = lm.lossmap_from_hits(hits)
hm = lm.CHitMap.from_hits(hits)
new_lm = hm.old_lossmap()
print(old_lm == new_lm)

ls = lm.CHitMap(settings.COLL_PATH)
lm.plot(ls)
# old_lm = lm.get_lossmap(settings.COLL_PATH)
# new_lm = ls.old_lossmap()
# print(old_lm == new_lm)

# ps = PhaseSpace(settings.STARTDIST_PATH)
# chs = ls.split(ps)
# lm.plot(chs)

# fig, ax = plt.subplots()
# ax.plot(ls.losses(), label='hits')
# ax.plot(ls.losses(True), label='integ')
# ax.legend(loc='upper right')
# # ax.set_xscale('log')
# plt.show()
