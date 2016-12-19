import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from collections import namedtuple
from math import pi

C_DIV_FREQ = 299793458.0/398765412.66*1e3 # converting to mm 

def phi_to_z(phi):
    return (pi - phi)/(2.0*pi)*C_DIV_FREQ

def comp():
    turns = []
    e_ramp = []

    BucketData = namedtuple("BucketData", 'e de phi z')
    six_part = BucketData([], [], [], [])
    toymodel = BucketData([], [], [], [])

    with open("1p.txt", 'r') as f:
        z_offset, de_offset = map(float, f.readline().strip().split())
        for line in f.readlines():
            if line.lstrip().startswith("#"): continue
            d = line.rstrip().split()
            turn = int(d[0])
            e_ref, e_part, z = map(float, d[1:])

            z -= z_offset
            e_part -= 0.45e6*de_offset # Kyrre suggested 0.5e6 here. Why the 0.5? Could it be because that is the diameter of the bucket?

            turns.append(turn)
            e_ramp.append(e_ref)
            six_part.e.append(e_part)
            six_part.de.append(e_part - e_ref)
            six_part.z.append(z)

    with open('toymodel_track.dat', 'r') as f:
        f.readline() # throw away first line

        turn_it = 0
        for turn in turns:
            energy = phase = 0
            while turn_it < turn:
                line = f.readline()
                energy, phase = map(float, line.rstrip().split(','))
                turn_it += 1
            toymodel.e.append(energy/1e6 + e_ramp[turn_it - 1])
            toymodel.de.append(energy/1e6)
            toymodel.z.append(phi_to_z(phase))


    # Energy plot
    fig, ax = plt.subplots()
    ax.plot(turns, six_part.e, color='b', zorder=1, label='sixtrack particle')
    ax.plot(turns, toymodel.e, color='red', zorder=2, label='toymodel particle')
    ax.plot(turns, e_ramp, color='black', zorder=10, label='reference energy')

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e3)))

    ax.legend(loc='lower right')
    ax.set_ylabel("Energy (GeV)")
    ax.set_xlabel("Turn")

    nbr_p = len(e_ramp)
    Point = namedtuple('Point', 'z e')
    six_center = Point(1.0/nbr_p*sum(six_part.z), 1.0/nbr_p*sum(six_part.de))
    toy_center = Point(1.0/nbr_p*sum(toymodel.z), 1.0/nbr_p*sum(toymodel.de))
    print("six center: ", six_center, "\ntoy center: ", toy_center)

    #Bucket diagram
    fig2, ax2 = plt.subplots()
    ax2.scatter(six_part.z, six_part.de, color='b', label='sixtrack')
    ax2.scatter(*six_center, color='b', marker='x', zorder=10)
    ax2.scatter(toymodel.z, toymodel.de, color='r', label='toymodel')
    ax2.scatter(*toy_center, color='r', marker='x', zorder=10)

    ax2.set_axisbelow(True)
    ax2.yaxis.grid(color='gray', linestyle='--')
    ax2.xaxis.grid(color='gray', linestyle='--')

    ax2.legend(loc='lower right')
    ax2.set_ylabel("∆E (MeV)")
    ax2.set_xlabel("z (mm)")

    fig2.suptitle("Bucket plot")

    plt.show()


comp()
