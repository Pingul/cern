import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def comp():
    turns = []
    e_part = []
    e_ramp = []

    toymodel = []

    with open("1p.txt", 'r') as f:
        for line in f.readlines():
            if line.lstrip().startswith("#"): continue
            turn, e_r, e_p = line.rstrip().split()
            turns.append(int(turn))
            e_ramp.append(float(e_r))
            e_part.append(float(e_p))

    with open('toymodel_track.dat', 'r') as f:
        f.readline() # throw away first line

        turn_it = 0
        for turn in turns:
            energy = 0
            while turn_it < turn:
                line = f.readline()
                energy = float(line.rstrip().split(',')[0])
                turn_it += 1
            toymodel.append(energy/1e6)


    e_toym = [e_ramp[i] + toymodel[i] for i, turn in enumerate(turns)]

    fig, ax = plt.subplots()
    ax.plot(turns, e_part, color='b', zorder=1, label='sixtrack particle')
    ax.plot(turns, e_toym, color='red', zorder=2, label='toymodel particle')
    ax.plot(turns, e_ramp, color='black', zorder=10, label='reference energy')

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e3)))

    ax.legend(loc='lower right')
    ax.set_ylabel("Energy (GeV)")
    ax.set_xlabel("Turn")
    plt.show()

comp()
