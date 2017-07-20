########
########
# This program compares a collimator motor function from a fill (received from Timber) with the corresponding trim function (from file). We are mostly
# concerned with the threshold interlocks. See 'coll_funcs.py' for what variables are used.
#
# Note: Will try to cache data in a directory 'cache' -- please ensure it exists and is empty, so no data is lost.
# 
# Use as:
# python3.5 verify_cmotors.py <action> (<argument>)
#
# Actions:
# - 'stats': Prints all variables for each collimator
# - 'hist': Creates a histogram over all collimators for a given variable.
#       <argument>: Variable to use, e.g. 'left_downstream'
# - 'plot': Plot motors and interlocks for a given collimator.
#       <arguemnt>: Collimator, e.g. 'TCP.D6L7.B1'
#
# Author: Joel Wretborn, CERN 12-07-2017
########

# SET THESE
# --------
# Timber fill used for the analysis
fill_nbr = 5699 

# Location from Trim variables
coll_trim_file = "motor_functions_ramp_B12_2017_CRS_2017_5_19_16_11_16.755108.csv" 

# Threshold in mm for allowed discrepancy between Timber/Trim
threshold = 0.1 

# Remove collimators from analysis
ignore_collimators = [ 
        "TCL.5R1.B1",
        "TCL.6R1.B1",
        "TCL.4R5.B1",
        "TCL.5R5.B1",
        "TCTPH.4L1.B1",
        "TCTPV.4L1.B1",
        "TCTPH.4L2.B1",
        "TCTPH.4L5.B1",
        "TCTPV.4L5.B1",
        "TCTPH.4L8.B1",
        "TCTPV.4L2.B1",
        "TCTPV.4L8.B1",
        "TCL.4L1.B2",
        "TCL.5L1.B2",
        "TCL.4L5.B2",
        "TCL.5L5.B2",
        "TCTPH.4R8.B2",
        "TCTPV.4R8.B2",
        "TCTPH.4R5.B2",
        "TCTPV.4R5.B2",
        "TCTPH.4R2.B2",
        "TCTPV.4R2.B2",
        "TCTPH.4R1.B2",
        "TCTPV.4R1.B2",
        "TCDQA.A4R6.B1",
        "TCDQA.A4L6.B2",
        ]
# --------

import sys
import numpy as np
import matplotlib.pyplot as plt
import coll_funcs as cf

def color_wrap(string, color):
    return '\x1b[{}m'.format(color) + string + '\x1b[0m'

def red(string):
    return color_wrap(string, "1;31;40")

def plot_func(cfunc):
    series = [
        ["t",
         "left_downstream",
         "dump_inner_left_downstream",
         "dump_outer_left_downstream"],
        ["t",
         "right_downstream",
         "dump_inner_right_downstream",
         "dump_outer_right_downstream"],
        ["t",
         "left_upstream",
         "dump_inner_left_upstream",
         "dump_outer_left_upstream"],
        ["t",
         "right_upstream",
         "dump_inner_right_upstream",
         "dump_outer_right_upstream"]]
    labels = ["LD", "RD", "LU", "RU"]
    
    fig, axs = plt.subplots(nrows=2, ncols=2)
    for tv, l, ax in zip(series, labels, axs.flatten()):
        d = cfunc.trim.get_series(names=tv)
        ax.plot(d[0], d[1], label='trim motor {}'.format(l), linestyle='--', zorder=4)
        ax.plot(d[0], d[2], label='trim int. in {}'.format(l), linestyle='--', zorder=4)
        ax.plot(d[0], d[3], label='trim int. out {}'.format(l), linestyle='--', zorder=4)

        df = cfunc.fill.get_series(names=cf.CollTrimVMap.fillVarsForTrim(tv))
        ax.plot(df[0], df[1], label='coll motor {}'.format(l))
        ax.plot(df[0], df[2], label='coll int. in {}'.format(l))
        ax.plot(df[0], df[3], label='coll int. out {}'.format(l))

        ax.legend(loc="upper right")
        ax.set_xlabel("t [s]")
        ax.set_ylabel("Pos. [mm]")
        ax.set_title("Collimator {} function {}".format(cfunc.fill.cname, l))

    # different plot for gaps
    fig2, ax = plt.subplots()
    tv = ["t", 
          "dump_inner_gap_upstream",
          "dump_inner_gap_downstream",
          "dump_outer_gap_upstream",
          "dump_outer_gap_downstream"]
    d = cfunc.trim.get_series(names=tv)
    df = cfunc.fill.get_series(names=cf.CollTrimVMap.fillVarsForTrim(tv))

    ax.plot(d[0], d[1], label='trim in GU'.format(l), linestyle='--', zorder=4)
    ax.plot(d[0], d[2], label='trim in GD'.format(l), linestyle='--', zorder=4)
    ax.plot(d[0], d[3], label='trim out GU'.format(l), linestyle='--', zorder=4)
    ax.plot(d[0], d[3], label='trim out GD'.format(l), linestyle='--', zorder=4)

    ax.plot(df[0], df[1], label='coll in GU'.format(l))
    ax.plot(df[0], df[2], label='coll in GD'.format(l))
    ax.plot(df[0], df[3], label='coll out GU'.format(l))
    ax.plot(df[0], df[3], label='coll out GD'.format(l))

    ax.legend(loc="upper right")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("Pos. [mm]")
    ax.set_title("Interlock gap for coll. {}".format(cfunc.fill.cname))

    plt.show()

def compare_funcs(funcs, threshold=0.1):
    tot_var = 0
    nbr_variables_exceeds = 0
    fills_exceeding = []

    for f in funcs:
        print(f.trim.cname)
        tv = f.trim.dmap[np.invert(f.trim.dmap == 't')]
        fv = cf.CollTrimVMap.fillVarsForTrim(tv)
        tv = cf.CollTrimVMap.trimVarsForFill(fv) # making sure everything is in the correct order

        td = f.trim.get_series(tv)
        fd = f.fill.get_series(fv)
        
        comp = np.amax(np.abs(td - fd), axis=1)

        fill_ok = True
        for v in zip(tv, fv, comp):
            tot_var += 1
            s = "\t{} - {} : max ∆={}".format(*v)
            if (v[2] > threshold):
                print(red(s))
                fill_ok = False
                nbr_variables_exceeds += 1
            else:
                print(s)

        if not fill_ok:
            fills_exceeding.append(f.fill.cname)
    
    print("\nStats:")
    print("\tThreshold: {} mm".format(threshold))
    print("\tVaribles exceeding threshold:       {}/{} ({:.0f}%)".format(nbr_variables_exceeds, tot_var, 100*float(nbr_variables_exceeds)/float(tot_var)))
    print("\tCollimators exceeding threshold:    {}/{} ({:.0f}%)".format(len(fills_exceeding), len(funcs), 100*float(len(fills_exceeding)/float(len(funcs)))))
    print("List of all collimators exceeding threshold")
    print("\t", end="")
    print("\n\t".join(fills_exceeding))
    
def plot_compare_hist(funcs, trim_var):
    delta_cg = np.empty(len(funcs))
    tv = np.array([trim_var])
    for i, f in enumerate(funcs):
        fv = cf.CollTrimVMap.fillVarsForTrim(tv)

        td = f.trim.get_series(tv)
        fd = f.fill.get_series(fv)
        comp = np.amax(np.abs(td - fd), axis=1)
        delta_cg[i] = comp[0]
    fig, ax = plt.subplots()
    bins = np.linspace(0, delta_cg.max(), 20)
    ax.hist(delta_cg, edgecolor='white', bins=bins)
    ax.set_xlabel("∆x [mm]")
    ax.set_ylabel("#")
    plt.title("Histogram for variable '{}'".format(trim_var))
    plt.show()


if __name__ == "__main__":
    action = sys.argv[1]
    if action == "dev":
        cf.cfunc_from_trimfile(coll_trim_file)
    else:

        funcs = cf.get_cfunctions(coll_trim_file, fill_nbr, ignore_collimators)
        print("\ninterpolate trims to match fills...", end=" ")
        for f in funcs:
            f.trim.recalculate_time(f.fill)
        print("done")

        if action == "stats":
            compare_funcs(funcs, threshold)
        elif action == "hist":
            variable = sys.argv[2]
            plot_compare_hist(funcs, variable)
        elif action == "plot":
            name = sys.argv[2]
            for f in funcs:
                if f.trim.cname == name:
                    plot_func(f)
                    break
            else:
                print("could not find collimator")
        else: 
            print("unrecognised action")
