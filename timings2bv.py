"""timings2bv.py

Create BrainVoyager .prt and .sdm files from a set of timing files.

One timing file represents one experimental condition and includes a single
line with times in msec (float), separated by spaces (e.g. 2.0 4.0 6.0).

Each timing file has the following naming convention:

    *_<condition_name>.* (e.g. myexp_condition1.txt)

"""


__author__ = 'Florian Krause <florian.krause@fladd.de>'
__version__ = '0.0.2'


import os
import random
import colorsys
from glob import glob

import numpy
import scipy.stats as stats


def _create_colours(n):
    colours = []
    for i in range(0, 360, 360/n):
        h = i / 360.0
        l = (50 + random.random() * 10) / 100
        s = (90 + random.random() * 10) / 100
        colours.append([int(x * 255) for x in colorsys.hls_to_rgb(h, l, s)])
    return colours

def hrf(tr, p=[6,16,1,1,6,0,32]):
    """Create an HRF from two gamma functions.

    Paramters
    ---------
    tr : float
        repetition time (in seconds) at which to sample the HRF

    p : list, optional
        parameters of the two gamma functions:
                                                            defaults
                                                            (seconds)
        p[0] - delay of response (relative to onset)         6
        p[1] - delay of undershoot (relative to onset)      16
        p[2] - dispersion of response                        1
        p[3] - dispersion of undershoot                      1
        p[4] - ratio of response to undershoot               6
        p[5] - onset (seconds)                               0
        p[6] - length of kernel (seconds)                   32


    """

    p = [float(x) for x in p]
    tr=float(tr)
    fMRI_T = 16.0

    # HRF in seconds
    dt = tr/fMRI_T
    u = numpy.arange(p[6] / dt + 1) - p[5] / dt
    g1 = stats.gamma.pdf(u, p[0] / p[2], scale=1.0 / (dt / p[2]))
    g2 = stats.gamma.pdf(u, p[1] / p[3], scale=1.0 / (dt / p[3])) / p[4]
    hrf = g1 - g2

    # Sample in volumes
    good_pts = numpy.array(range(numpy.int((p[6] + tr) / tr))) * fMRI_T
    hrf = hrf[[int(x) for x in list(good_pts)]]
    hrf = hrf / numpy.sum(hrf);

    return hrf

def convert2prt(input_mask, run_length, tr, msec=False, durations=None, output=None):
    if durations is None:
        if not msec:
            durations = 0
        else:
            durations = 0.001
    if output is None:
        output = "Simulation"

    files = glob(os.path.join(input_mask))
    with open(output + ".prt", 'w') as f:
        f.write("\n")
        f.write("FileVersion:        3\n")
        f.write("\n")
        if not msec:
            f.write("ResolutionOfTime:   Volume\n")
        else:
            f.write("ResolutionOfTime:   msec\n")
        f.write("\n")
        f.write("Experiment:         {0}\n".format(output))
        f.write("\n")
        f.write("BackgroundColor:    0 0 0\n")
        f.write("TextColor:          255 255 255\n")
        f.write("TimeCourseColor:    255 255 255\n")
        f.write("TimeCourseThick:    4\n")
        f.write("ReferenceFuncColor: 0 0 80\n")
        f.write("ReferenceFuncThick: 3\n")
        f.write("\n")
        f.write("ParametricWeights:  0\n")
        f.write("\n")
        f.write("NrOfConditions:     {0}\n".format(len(files)))
        if not msec:
            rjust = 4
        else:
            rjust = 8
        colours = _create_colours(len(files))
        for cnt, condition in enumerate(files):
            with open(condition) as c:
                data = c.readline()
                events = [float(x) for x in data.strip().split(" ")]
                f.write("\n")
                f.write(condition.split("_")[-1].split(".")[0] + "\n")
                f.write(repr(len(events)) + "\n")
                for e, event in enumerate(events):
                    if not msec:
                        begin = int(event / tr) + 1
                        end = int(begin + (durations / tr) - 1)
                    else:
                        begin = int(event * 1000)
                        end = int(begin + durations * 1000)
                    f.write("{0} {1}\n".format(
                        repr(begin).rjust(rjust),
                        repr(end).rjust(rjust)))
                f.write("Color: {0} {1} {2}\n".format(colours[cnt][0],
                                                      colours[cnt][1],
                                                      colours[cnt][2]))

def convert2sdm(input_mask, run_length, tr, msec=False, durations=None, output=None):
    if durations is None:
        if msec:
            durations = 0.001
        else:
            durations = 0
    if output is None:
        output = "Simulation"

    files = glob(os.path.join(input_mask))
    regressors = []
    for file in files:
        with open(file) as f:
            data = f.readline()
            regressors.append([float(x) for x in data.strip().split(" ")])

    if not msec:
        convolved = []
        for regressor in regressors:
            tmp = [0] * int(run_length / tr)
            for c,x in enumerate(regressor):
                idx = int(x / tr)
                for d in range(idx, idx + int(durations / tr)):
                    tmp[d] = 1
            hrf_ = hrf(tr)
            cregressor = numpy.convolve(tmp, hrf_)[0:len(tmp)]
            cregressor = ["{:.6f}".format(round(x, 6)) for x in cregressor]
            convolved.append(cregressor)

    else:
        convolved = []
        for regressor in regressors:
            tmp = [0] * (run_length * 1000)
            for c,x in enumerate(regressor):
                idx = int(float(x) * 1000)
                for d in range(idx, idx + int(durations * 1000)):
                    tmp[d] = 1
            hrf_ = hrf(0.001)
            cregressor = numpy.convolve(tmp, hrf_)[0:len(tmp)]
            cregressor = cregressor[0::int(tr*1000)]
            cregressor = ["{:.6f}".format(round(x, 6)) for x in cregressor]
            convolved.append(cregressor)

    with open(output + ".sdm", 'w') as f:
        f.write("FileVersion:            1\n")
        f.write("\n")
        f.write("NrOfPredictors:         {0}\n".format(len(files) + 1))
        f.write("NrOfDataPoints:         {0}\n".format(run_length / tr))
        f.write("IncludesConstant:       1\n")
        f.write("FirstConfoundPredictor: {0}\n".format(len(files) + 1))
        f.write("\n")
        colours = _create_colours(len(files))
        f.write("   ".join(
            [str(x).strip("[]").replace(",","") for x in colours]) + \
                "   255 255 255\n")
        names = [x.split("_")[-1].split(".")[0] for x in files]
        f.write(" ".join(['"{0}"'.format(x) for x in names]) + ' "Constant"\n')
        for c,_ in enumerate(convolved[0]):
            f.write(" ".join(
                [regressor[c].rjust(10) for regressor in convolved]) + \
                    "   1.000000\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert timing files to .prt/.sdm")
    parser.add_argument('input_mask', type=str,
                        help="mask of input files (e.g. 'stimes/stimes.001_*')")
    parser.add_argument('run_length', type=int,
                        help='run length in seconds')
    parser.add_argument('tr', type=float,
                        help='repetition time in seconds')
    parser.add_argument('-m', '--msec', action="store_true",
                        help='use msec-based timing')
    parser.add_argument('-d', '--durations', type=float,
                        help='durations of events in seconds (if epochs)')
    parser.add_argument('-o', '--output', type=str,
                        help='name of resulting output files')
    args = parser.parse_args()

    convert2prt(args.input_mask, args.run_length, args.tr, args.msec,
                args.durations, args.output)
    convert2sdm(args.input_mask, args.run_length, args.tr, args.msec,
                args.durations, args.output)
