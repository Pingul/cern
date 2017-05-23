import sys
sys.path.insert(0, "../common")
import numpy as np
import matplotlib.pyplot as plt
import coll_funcs as cf

def color_wrap(string, color):
    return '\x1b[{}m'.format(color) + string + '\x1b[0m'

def red(string):
    return color_wrap(string, "1;31;40")

fill_nbr = 5668
coll_trim_file = "motor_functions_ramp_B12_2017_CRS_2017_5_19_16_11_16.755108.csv"
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

def plot_func(cfunc):
    trim_vars = [ 
        ["t",
        "left_downstream",
        "right_downstream",
        "left_upstream",
        "right_upstream",],
        ["t",
        "dump_inner_left_downstream",
        "dump_inner_left_upstream",
        "dump_inner_right_downstream",
        "dump_inner_right_upstream",],
        ["t",
        "dump_outer_left_downstream",
        "dump_outer_left_upstream",
        "dump_outer_right_downstream",
        "dump_outer_right_upstream",]]

    for tv in trim_vars:
        fig, ax = plt.subplots(nrows=2, ncols=1)
        
        d = cfunc.trim.get_series(names=("t", "left_downstream", "right_downstream", "left_upstream", "right_upstream"))
        ax[0].plot(d[0], d[1], label='trim LD', linestyle='--', zorder=4)
        ax[0].plot(d[0], d[2], label='trim RD', linestyle='--', zorder=4)
        ax[1].plot(d[0], d[3], label='trim LU', linestyle='--', zorder=4)
        ax[1].plot(d[0], d[4], label='trim RU', linestyle='--', zorder=4)
        
        df = cfunc.fill.get_series(names=("t", "MEAS_MOTOR_LD", "MEAS_MOTOR_RD", "MEAS_MOTOR_LU", "MEAS_MOTOR_RU"))
        ax[0].plot(df[0], df[1], label='fill LD')
        ax[0].plot(df[0], df[2], label='fill RD')
        ax[1].plot(df[0], df[3], label='fill LU')
        ax[1].plot(df[0], df[4], label='fill RU')

        for a in ax:
            a.legend(loc="upper right")
            a.set_xlabel("t [s]")
            a.set_ylabel("Collimator pos. [mm]")
        
        fig.suptitle("Collimator function {}".format(cfunc.fill.cname))
    plt.show()

def compare_funcs(funcs, threshold=0.1):
    print("Compare")

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
    plt.title("Compare collgaps for '{}'".format(trim_var))
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
            compare_funcs(funcs, 0.1)
        elif action == "hist":
            variable = sys.argv[2]
            plot_compare_hist(funcs, variable)
        elif action == "plot-fill":
            plot_func(funcs[-1])
        else: 
            print("unrecognised action")
