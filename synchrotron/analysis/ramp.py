import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

def read_ramp(filePath, nbr_turns):
    e = np.empty(nbr_turns)
    with open(filePath, 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i >= nbr_turns: break
            e[i] = float(line.rstrip().split()[1])*1e6
    return e

    # turns = range(nbr_turns)
    # #de = np.gradient(e)
    # de = np.diff(e) # size n - 1
    # de = np.append(de, de[-1])

    # slope, intercept, r_value, p_value, std_err = stats.linregress(turns, de)
    # de_fitted = [slope*turn + intercept for turn in turns]
    # e_fitted = []
    # s = 0
    # for v in de_fitted:
        # s += v
        # e_fitted.append(s + 450e9)

    # return {'e' : e, 'e_fitted' : e_fitted, 'de' : de, 'de_fitted' : de_fitted}

def plot_ramps():
    ramp_files = [ "10s_linear_ramp.dat", "20s_linear_ramp.dat", "30s_linear_ramp.dat", "40s_linear_ramp.dat", "50s_linear_ramp.dat", "interpolated_ramp.dat", "LHC_ramp.dat"]

    fig, ax = plt.subplots()
    nbr_turns = 50*11245
    turns = range(nbr_turns)

    color_list = plt.cm.Set3(np.linspace(0, 1, len(ramp_files)))
    for f in ramp_files:
        ax.plot(turns, read_ramp("resources/{}".format(f), nbr_turns), label=f)

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.2f}".format(x/11245.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.0f}".format(x/1e9)))
    ax.set_ylabel("Energy (GeV)")
    ax.set_xlabel("t (s)")
    ax.legend(loc="upper right")
    plt.title("Ramping functions")
    plt.show()


if __name__ == "__main__":
    plot_ramps()


