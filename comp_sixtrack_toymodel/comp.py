import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from collections import namedtuple
from math import pi

C_DIV_FREQ = 299793458.0/400e6*1e3 # converting to mm 

def phi_to_z(phi):
    return (phi - pi)/(2.0*pi)*C_DIV_FREQ

def comp():
    turns = []
    e_ramp = []

    BucketData = namedtuple("BucketData", 'e de phi')
    six_part = BucketData([], [], [])
    toymodel = BucketData([], [], [])

    with open("1p.txt", 'r') as f:
        for line in f.readlines():
            if line.lstrip().startswith("#"): continue
            d = line.rstrip().split()
            turn = int(d[0])
            e_ref, e_part, e_phi = map(float, d[1:])

            turns.append(turn)
            e_ramp.append(e_ref)
            six_part.e.append(e_part)
            six_part.de.append(e_part - e_ref)
            six_part.phi.append(e_phi)

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
            toymodel.phi.append(phi_to_z(phase))


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
    Point = namedtuple('Point', 'x y')
    six_center = Point(1.0/nbr_p*sum(six_part.phi), 1.0/nbr_p*sum(six_part.de))
    toy_center = Point(1.0/nbr_p*sum(toymodel.phi), 1.0/nbr_p*sum(toymodel.de))
    print("six center: ", six_center, "toy center: ", toy_center)

    #Bucket diagram
    fig2, ax2 = plt.subplots()
    ax2.scatter(six_part.phi, six_part.de, color='b', label='sixtrack')
    ax2.scatter(*six_center, color='b', marker='x')
    ax2.scatter(toymodel.phi, toymodel.de, color='r', label='toymodel')
    ax2.scatter(*toy_center, color='r', marker='x')

    ax2.set_axisbelow(True)
    ax2.yaxis.grid(color='gray', linestyle='--')
    ax2.xaxis.grid(color='gray', linestyle='--')

    ax2.legend(loc='lower right')
    ax2.set_ylabel("∆E (MeV)")
    ax2.set_xlabel("z (mm)")

    fig2.suptitle("Bucket plot")

    plt.show()


comp()
