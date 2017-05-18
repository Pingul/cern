import numpy as np
from oml import Fill, fills_from_file

import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from scipy import stats, interpolate
import math

def moving_average(sequence, N):
    """ Moving average given the sequence. Returns an list equal in length
    to the one given """

    average = np.convolve(sequence, np.ones((N,))/N, mode='same')
    return average

def merge_variable_interpolated(fills, var, resolution=0.1):
    xmin = -sys.maxsize
    xmax = -xmin
    for fill in fills:
        xmin = max(fill.data[var].x.min(), xmin)
        xmax = min(fill.data[var].x.max(), xmax)

    newx = np.arange(math.ceil(xmin), math.floor(xmax), resolution)

    avg = np.zeros(newx.size)
    maxy = np.zeros(newx.size)
    miny = np.full(maxy.shape, sys.maxsize, dtype=float)

    for fill in fills:
        ip = interpolate.interp1d(*fill.data[var], assume_sorted=True)(newx)
        avg += ip
        maxy = np.maximum(maxy, ip)
        miny = np.minimum(miny, ip)
    avg /= len(fills)
    return {
            "average" : Fill.Variable((newx, avg)),
            "max" : Fill.Variable((newx, maxy)),
            "min" : Fill.Variable((newx, miny))
            }


def merge_variable(fills, var):
    return merge_variable_interpolated(fills, var)
    xmin = sys.maxsize
    xmax = -xmin
    for fill in fills:
        xmin = min(fill.data[var].x.min(), xmin)
        xmax = max(fill.data[var].x.max(), xmax)

    # The x-series
    newx = np.arange(xmin, xmax + 1, 1, dtype=np.int) # using int's to maybe make comparisons more safe

    # Will store the average values of y
    avgy = np.zeros(newx.size, dtype=np.float)
    norm = np.zeros(newx.size, dtype=np.float) # we use this to basically take the average value at each timestamp

    # Will store maximum and minimum values of y
    maxy = np.zeros(newx.size, dtype=np.float)
    miny = 1.0e10*np.ones(newx.size, dtype=np.float)

    for fill in fills:
        x_i = fill.data[var].x.astype(np.int)
        index_map = np.searchsorted(newx, x_i)

        np.add.at(avgy, index_map, fill.data[var].y)
        np.add.at(norm, index_map, 1.0)

        np.maximum.at(maxy, index_map, fill.data[var].y)
        np.minimum.at(miny, index_map, fill.data[var].y)

    # delete empty values (where no data was written)
    to_delete = np.where(norm < 0.5)
    newx = np.delete(newx, to_delete)
    maxy = np.delete(maxy, to_delete)
    miny = np.delete(miny, to_delete)
    norm = np.delete(norm, to_delete)
    avgy = np.delete(avgy, to_delete)
    avgy /= norm

    newx = newx.astype(dtype=float)
    return {
            "average" : Fill.Variable((newx, avgy)),
            "max" : Fill.Variable((newx, maxy)),
            "min" : Fill.Variable((newx, miny))
            }


def aggregate_fill(beam, fill_list=[], from_cache=False):
    """ Create an aggregate fill of the fill ids. If reading from cache, 
        the fill_list can be empty. """

    if not fill_list and not from_cache:
        raise Exception("'fill_list' can't be empty if not to read from cache'")
    elif not beam in (1, 2):
        raise Exception("'beam' must be 1 or 2, is '{}'".format(beam))

    fill = Fill(nbr=-beam, beam=beam, fetch=False)
    if from_cache:
        fill.load_cache()
        return fill

    fills = []
    for nbr in fill_list:
        _fill = Fill(nbr)
        _fill.beta_coll_merge()
        fills.append(_fill)
        
    # This is cumbersome to refactor... TODO
    var = ["synch_coll_b{}".format(beam), "beta_coll_b{}".format(beam), "energy", "intensity_b{}".format(beam)]
    for v in var:
        merged = merge_variable(fills, v)
        fill.data[v] = merged["average"]
    return fill

def plot_aggregate_fill_overlayed(beam, fill_list):
    fig, ax = plt.subplots()
    for nbr in fill_list:
        fill = Fill(nbr, beam=beam)
        ax.plot(*fill.blm_ir3(), color=np.random.rand(3), alpha=0.3)
        # ax.plot(*fill.energy(), color=np.random.rand(3), alpha=0.3)
        # ax.plot(fill.motor_ir7().x, fill.motor_ir7().y[1], color=np.random.rand(3), alpha=0.3)

    aggr = aggregate_fill(beam, fill_list)
    ax.plot(*aggr.blm_ir3(), color='black', label='Aggregate', zorder=5)
    ax.set_xlim(aggr.blm_ir3().x[aggr.OML_period()] + np.array([-5, +120]))
    ax.set_yscale("log")

    ax.set_xlabel("t (s)")
    ax.set_ylabel("TCP IR3 BLM signal")
    plt.title("Overlay plot (beam {})".format(beam))
    plt.show()



def plot_aggregate_fill(beam, fill_list):
    # We can't reuse the 'aggregate_fill' function above, as we want to plot more
    # than for a normal fill (min, max, average)
    fills = []
    for nbr in fill_list:
        fill = Fill(nbr, beam=beam)
        fill.beta_coll_merge()
        fills.append(fill)
        
    var = ["synch_coll_b{}".format(beam), "beta_coll_b{}".format(beam), "energy"]
    merged = {}
    for v in var:
        merged[v] = merge_variable(fills, v)

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    e_ax = ax.twinx()

    v = "synch_coll_b{}".format(beam)
    ax.plot(*merged[v]["average"], zorder=10, color='red', label='OML')
    ax.fill_between(*merged[v]["max"], merged[v]["min"].y, color='red', alpha='0.3')

    v = "beta_coll_b{}".format(beam)
    ax.plot(*merged[v]["average"], zorder=9, color='green', label='TL')
    ax.fill_between(*merged[v]["max"], merged[v]["min"].y, color='green', alpha='0.2')

    e_ax.plot(*merged["energy"]["average"], zorder=5, color='black', label='energy')
    e_ax.set_ylabel("Reference energy (GeV)")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("BLM signal")
    ax.legend(loc="upper right")
    ax.set_xlim(fills[0].blm_ir3().x[fills[0].OML_period()] + np.array([-5, +120]))
    # ax.axvspan(0.0, 13.55, facecolor='b', zorder=0, alpha=0.1)
    plt.title("Aggregate fill 2016 (beam {})".format(beam))
    plt.show()

#############
## Statistics

def draw_histogram(title, data, binsize, xlabel='', ylabel=''):
    maxbin = max(data) + binsize
    minbin = min(data)
    bins = np.arange(minbin, maxbin, binsize)
    fig, ax = plt.subplots()
    ax.hist(data, bins=bins, edgecolor='white')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()


def histogram_intensity(file):
    fills = fills_from_file(file, "OML")
    intensities = []
    for nbr in fills:
        fill = Fill(nbr, fetch=False)
        fill.fetch()
        intensities.append(max(fill.intensity().y))

    draw_histogram('Intensity for {}'.format(file), intensities, 1e13, "Intensity", "Count")

def histogram_spike_duration(file):
    fills = fills_from_file(file, "OML")
    outliers = []
    durations = []
    for nbr in fills:
        fill = Fill(nbr)

        start, end = fill.OML_period()
        d = fill.blm_ir7().x[end] - fill.blm_ir7().x[start]
        if d < 70 or d > 300:
            outliers.append(nbr)
        durations.append(d)

    draw_histogram('Spike duration for {}'.format(file), durations, 10, 'Seconds', 'Count')
    return outliers

def histogram_crossover_point(file):
    """ Histogram on what durations OML dominate the transversal losses """
    fills = fills_from_file(file, "OML")
    durations = []
    outliers = []
    for nbr in fills:
        fill = Fill(nbr)
        
        crossover = fill.crossover_point()
        if crossover['t'] > 40:
            outliers.append(nbr)
        else:
            durations.append(crossover["t"])

    draw_histogram("Duration OML > transversal losses from '{}'".format(file), durations, 1, 'Duration (s) after spike with OML > TM', 'Count')
    return outliers

def histogram_max_spike(fills):
    spike_time = {1:[], 2:[]}

    for nbr in fills:
        for b in (1, 2):
            fill = Fill(nbr, beam=b)
            losses = fill.blm_ir3()
            ispike = np.argmax(losses.y)
            spike_time[b].append(losses.x[ispike])
    fig, ax = plt.subplots()
    ax.hist([spike_time[1], spike_time[2]], label=["Beam 1", "Beam 2"])
    ax.legend(loc="upper right")
    ax.set_xlabel("Delta t (s) from start of ramp till maximum OML")
    ax.set_ylabel("Fill count")
    plt.show()
    # draw_histogram('Max spike event for beam {}'.format(beam), spike_time, 0.5, 'Delta t (s) from start of ramp till spike', 'Count')

def histogram_max_abort_gap_before_OML(file):
    fills = fills_from_file(file, "OML")
    max_ag = {1:[], 2:[]}

    for nbr in fills:
        for beam in (1, 2):
            fill = Fill(nbr, beam=beam)
            OML_start = fill.OML_period()[0]
            max_ag[beam].append(fill.abort_gap().y[:OML_start].max())
    bins = np.arange(0, 2e10, 2e9)
    fig, ax = plt.subplots()
    ax.hist([max_ag[1], max_ag[2]], bins=bins, label=["Beam 1", "Beam 2"])
    ax.legend(loc="upper right")
    ax.set_xlabel("Max abort gap intensity before start of ramp")
    ax.set_ylabel("Fill count")
    plt.show()

def histogram_OML_peak(file):
    fills = fills_from_file(file, "OML")
    OML_peak = {1: [], 2: []}

    for nbr in fills:
        for beam in (1, 2):
            fill = Fill(nbr, beam=beam)
            OML_peak[beam].append(np.log10(fill.blm_ir3().y.max()))
    fig, ax = plt.subplots()
    ax.hist([OML_peak[1], OML_peak[2]], label=["Beam 1", "Beam 2"])
    ax.legend(loc="upper right")
    ax.set_xlabel("BLM peak signal [$10^x$]")
    ax.set_ylabel("Fill count")
    plt.title("Peak IR3 TCP BLM signal distribution")
    plt.show()


def bar_graph_crossover_point(file):
    """ Bar graph displaying the 'time till max spike' + 'time where OML > TM' for all fills in file """
    fills = fills_from_file(file, "OML")
    oml_dom_duration = []
    time_till_spike = []
    total_duration = []

    odd = []
    for nbr in fills:
        fill = Fill(nbr)

        crossover = fill.crossover_point()
        oml = fill.blm_ir3()
        ispike = np.argmax(oml.y)
        tspike = oml.x[ispike]
        time_till_spike.append(tspike)

        oml_dom_duration.append(crossover["t"] - tspike)

        total_duration.append(crossover['t'])

        if crossover['t'] > 40 or tspike > 12: # why am I using this?
            odd.append(nbr)

    print("odd looking fills: ", odd)
    fig, ax = plt.subplots()
    ax.bar(range(len(oml_dom_duration)), time_till_spike, label='Time to max peak (s)')
    ax.bar(range(len(oml_dom_duration)), oml_dom_duration, bottom=time_till_spike, label='Crossover point (s)')
    ax.legend(loc="upper right")
    ax.set_xlabel("Fill nbr")
    ax.set_ylabel("Duration (s)")
    plt.show()

def comp_blm_ir3_vs_ir7(file):
    fills = fills_from_file(file, "OML")
    ok = 0
    notok = 0
    sdata = {
        'max' : [],
        'mean' : [],
    }
    bdata = {
        'max' : [],
        'mean' : []
    }
    for nbr in fills:
        fill = Fill(nbr)
        fill.beta_coll_merge()

        smin, smax = fill.OML_period()
        tmin, tmax = fill.blm_ir3().x[[smin, smax]]
        bmin = fill.blm_ir7().index_for_time(tmin)
        bmax = fill.blm_ir7().index_for_time(tmax)

        bsubset = fill.blm_ir7().y[bmin:bmax]
        ssubset = fill.blm_ir3().y[smin:smax]

        sdata['max'].append(max(ssubset))
        sdata['mean'].append(np.mean(ssubset))
        bdata['max'].append(max(bsubset))
        bdata['mean'].append(np.mean(bsubset))

    fig, ax = plt.subplots()
    ax.set_xlabel("Synchrotron (IR3) TCP")
    ax.set_ylabel("Betatron (IR7) TCPs")

    ax.scatter(sdata['max'], bdata['max'], color='r')
    slope, intercept, r_value, p_value, std_err = stats.linregress(sdata['max'], bdata['max'])
    # print(slope, intercept, r_value, p_value, std_err)
    xval = [0, 1]
    max_yval = [slope*x + intercept for x in xval]
    ax.plot(xval, max_yval, color='r', label='max')

    ax.scatter(sdata['mean'], bdata['mean'], color='b')
    slope, intercept, r_value, p_value, std_err = stats.linregress(sdata['mean'], bdata['mean'])
    # print(slope, intercept, r_value, p_value, std_err)
    mean_yval = [slope*x + intercept for x in xval]
    ax.plot(xval, mean_yval, color='b', label='mean')

    ax.plot([0, 1], [0, 1], color = 'black', label='delimiter')

    for v in ['max', 'mean']:
        count = 0
        for i, sd in enumerate(sdata[v]):
            if bdata[v][i] > sd: 
                count += 1
        print(v, "over: ", count, "({}%)".format(int(float(count)/len(sdata[v])*100)))

    plt.title('Losses due to synchrotron vs betatron oscillations\n for {}'.format(file))
    ax.legend(loc='upper right')
    ax.set_ylim([0, 0.5])
    ax.set_xlim([0, 0.5])
    plt.show()

def comp_blm_ir3_vs_intensity(file):
    fills = fills_from_file(file, "OML")
    intensity = []
    mean_loss = []
    max_loss = []
    discarded = 0
    for nbr in fills:
        fill = Fill(nbr, False)
        fill.fetch()
        smin, smax = fill.OML_period()
        ssubset = fill.blm_ir3().y[smin:smax]

        maxint = max(fill.intensity().y)
        if maxint < 1.8e14:
            discarded += 1
            continue

        mean_loss.append(np.mean(ssubset))
        max_loss.append(max(ssubset))
        intensity.append(maxint)

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharey=ax1) 

    ax1.set_xlabel("Mean momentum (IR3) TCP")
    ax1.set_ylabel("Intensity")
    ax1.scatter(mean_loss, intensity, color='b', label='mean')
    ax1.set_xlim([0, 1.1*max(mean_loss)])
    ax1.set_ylim([1.5e14, 1.1*max(intensity)])
    ax1.legend(loc="lower right")

    ax2.set_xlabel("Max momentum (IR3) TCP")
    ax2.set_ylabel("Intensity")
    ax2.scatter(max_loss, intensity, color='r', label='max')
    ax2.set_xlim([0, 1.1*max(max_loss)])
    ax2.legend(loc="lower right")

    percent_used = int(round(float(len(intensity))/(len(intensity) + discarded) * 100))
    fig.suptitle("Intensity vs OML for {} (only intenities > 1.8e14, {}% of total)\n".format(file, percent_used))

    plt.show()

def comp_blm_ir3_vs_abort_gap(file):
    fills = fills_from_file(file, "OML")
    abort_gap = []
    average_loss = []
    max_loss = []
    for nbr in fills:
        fill = Fill(nbr, False)
        fill.fetch()
        smin, smax = fill.OML_period()

        # Only looking until t_co instead -- will not affect max
        smax = fill.crossover_point()['i']

        tmax = fill.blm_ir3().x[smax]
        tmin = fill.blm_ir3().x[smin]

        # tmax = find_crossover_point(fill)['t']

        ag_average = moving_average(fill.abort_gap().y, 5)
        agmin = fill.abort_gap().index_for_time(tmin)
        agmax = fill.abort_gap().index_for_time(tmax)

        ssubset = fill.blm_ir3().y[smin:smax]

        average_loss.append(np.average(ssubset))
        max_loss.append(max(ssubset))
        abort_gap.append(ag_average[agmin] - ag_average[agmax])


    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharey=ax1) 

    # fig1, ax1 = plt.subplots()
    ax1.set_xlabel("Average BLM")
    ax1.set_ylabel("∆ abort gap intensity")
    ax1.scatter(average_loss, abort_gap, color='b', label='average')
    ax1.set_xlim([0, 1.1*max(average_loss)])
    ax1.set_ylim([0, 1.1*max(abort_gap)])

    xval = [0, 1]
    slope, intercept, r_value, p_value, std_err = stats.linregress(average_loss, abort_gap)
    print("Average fit")
    print("\tk  ={:>10.3E}\n\tm  ={:>10.3E}\n\tr  ={:>10.7f}\n\tp  ={:>10.3E}\n\te^2={:>10.3E}".format(slope, intercept, r_value, p_value, std_err))
    yfit = [slope*x + intercept for x in xval]
    ax1.plot(xval, yfit, color='gray')

    ax1.legend(loc="lower right")

    # fig2, ax2 = plt.subplots()
    ax2.set_xlabel("Max BLM")
    ax2.scatter(max_loss, abort_gap, color='r', label='max')
    ax2.set_xlim([0, 1.1*max(max_loss)])
    ax2.legend(loc="lower right")

    slope, intercept, r_value, p_value, std_err = stats.linregress(max_loss, abort_gap)
    print("Max fit")
    print("\tk  ={:>10.3E}\n\tm  ={:>10.3E}\n\tr  ={:>10.7f}\n\tp  ={:>10.3E}\n\te^2={:>10.3E}".format(slope, intercept, r_value, p_value, std_err))
    yfit = [slope*x + intercept for x in xval]
    ax2.plot(xval, yfit, color='gray')

    fig.suptitle("Correlation between abort gap intensity and BLM signal for TCP in IR3")
    plt.show()

def corr_inj_OML(fills):
    inj_t = {"b1" : [], "b2" : []}
    oml = {"b1" : [], "b2" : []}
    for i, nbr in enumerate(fills):
        for b in (1, 2):
            f = Fill(nbr, beam=b)

            inj = next((item for item in f.meta['beamModes'] if item['mode'] == 'INJPHYS'), None)
            preramp = next((item for item in f.meta['beamModes'] if item['mode'] == 'PRERAMP'), None)

            if preramp is None or inj is None:
                raise Exception("Invalid fill {}".format(nbr))

            if b == 1 and f.blm_ir3().y.max() > 0.3 or b == 2 and f.blm_ir3().y.max() > 0.11:
                continue


            dt = preramp['endTime'] - inj['startTime']
            inj_t['b{}'.format(b)].append(dt)
            oml['b{}'.format(b)].append(f.blm_ir3().y.max())

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, ax = plt.subplots(nrows=2, ncols=1)

    ax[0].scatter(inj_t['b1'], oml['b1'], label='beam 1')
    k1, m1, r1, p1, err1 = stats.linregress(inj_t['b1'], oml['b1'])
    ax[0].plot(ax[0].get_xlim(), [k1*x + m1 for x in ax[0].get_xlim()], label='fit b1, r={:.3f}'.format(r1))

    ax[1].scatter(inj_t['b2'], oml['b2'], label='beam 2', c=colors[1])
    k2, m2, r2, p2, err2 = stats.linregress(inj_t['b2'], oml['b2'])
    ax[1].plot(ax[1].get_xlim(), [k2*x + m2 for x in ax[1].get_xlim()], label='fit b2, r={:.3f}'.format(r2), c=colors[1])

    for a in ax:
        a.legend(loc="upper right")
        a.set_xlabel("t")
        a.set_ylabel("OML peak")
    print("Fit 1: {}*x + {}, r={}".format(k1, m1, r1))
    print("Fit 2: {}*x + {}, r={}".format(k2, m2, r2))
    plt.show()

