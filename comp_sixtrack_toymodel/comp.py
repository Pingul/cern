import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from collections import namedtuple
from math import pi, sqrt
import numpy as np

def phi_to_z(phi, tot_energy):
    """ tot_energy in MeV """
    C = 26658.8832
    h = 35640.0
    b = sqrt(1.0 - (938.2796/tot_energy)**2)
    z = (pi - phi)*C/(2*pi*h*b)
    return z*1e3 # convert to mm

# F_RF = 400788731.3727857
# C_DIV_FREQ = 299793458.0/F_RF*1e3 # converting to mm 
# def phi_to_z(phi):
    # return (pi - phi)/(2.0*pi)*C_DIV_FREQ

def imax(data):
    """ Returns both the max and the index for the max value """
    i = max(range(len(data)), key = data.__getitem__)
    return (data[i], i)

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
            e_part -= 0.45e6*de_offset 

            turns.append(turn)
            e_ramp.append(e_ref)
            six_part.e.append(e_part)
            six_part.de.append(e_part - e_ref)
            six_part.z.append(z)

    turns = np.array(turns)

    with open('toymodel_track.dat', 'r') as f:
        f.readline() # throw away first line

        turn_it = 0
        for turn in turns:
            energy = phase = 0
            while turn_it < turn:
                line = f.readline()
                if len(line) == 0: break
                energy, phase, h = map(float, line.rstrip().split(','))
                turn_it += 1
            toymodel.e.append(energy/1e6 + e_ramp[turn_it - 1])
            toymodel.de.append(energy/1e6)
            total_energy = toymodel.e[-1] + toymodel.de[-1]
            toymodel.z.append(phi_to_z(phase, total_energy))

    nbr_p = len(e_ramp)

    # Frequency
    N = nbr_p
    six_freq = np.abs(np.fft.fft(six_part.e)[:N])/N
    toy_freq = np.abs(np.fft.fft(toymodel.e)[:N])/N

    six_freq[0] = toy_freq[0] = 0.0 # Removing the 'constant' frequency
    msix, isix = imax(six_freq)
    mtoy, itoy = imax(toy_freq)

    f = np.fft.fftfreq(turns.shape[-1])
    print("Frequency\n\tsix=", f[isix], ", toy=", f[itoy])

    if isix != itoy:
        print("Frequencies diverge. Plotting frequencies...")

        fftfig, fftax = plt.subplots()
        fftax.plot(f, six_freq, c='b', label='sixtrack')
        fftax.plot(f, toy_freq, c='r', label='toymodel')

        fftax.legend(loc='upper right')
        fftax.set_xlabel("frequency (w.r.t turns)")
        fftfig.suptitle("Frequency of synchrotron oscillation")


    # Energy plot
    fig, ax = plt.subplots()
    ax.plot(turns, six_part.e, color='b', zorder=1, label='sixtrack particle')
    ax.plot(turns, toymodel.e, color='red', zorder=2, linestyle='--', label='toymodel particle')
    ax.plot(turns, e_ramp, color='black', zorder=10, label='reference energy')

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e3)))

    ax.legend(loc='lower right')
    ax.set_ylabel("Energy (GeV)")
    ax.set_xlabel("Turn")

    # voltage = [6 + 2.9491187074838457087e-07*t for t in turns]
    # vax = ax.twinx()
    # vax.get_yaxis().set_visible(False)
    # vax.plot(turns, voltage, color='gray', zorder=0, label='voltage')

    Point = namedtuple('Point', 'z e')
    six_center = Point(1.0/nbr_p*sum(six_part.z), 1.0/nbr_p*sum(six_part.de))
    toy_center = Point(1.0/nbr_p*sum(toymodel.z), 1.0/nbr_p*sum(toymodel.de))
    print("Center\n\tsix=", six_center, "\n\toy=", toy_center)

    # Bucket diagram
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


if __name__ == "__main__":
    comp()
