Small project aimed to help understand the longitudinal motion of particles inside a synchrotron accelerator. 

# Compiling

The code has been tested on macOS Sierra and on LXPLUS. To compile it the first time, simply write `make` in the main directory.

Options:
- It is possible to run the code multi-threaded. Go into the `makefile` and change `MULTI_THREADING=1`. Note that [Threading Building Blocks (TBB)](https://www.threadingbuildingblocks.org) is needed for multi-threading. Installing it on LXPLUS is a little bit tricky, so in general I recommend running the single threaded version there. Installing TBB on macOS is easy using [MacPorts](https://www.macports.org), and it should be possible to use `apt-get` on Linux.

# Use the code
There are 2 main use cases for the tool:

### 1. Simulation meant to be analysed using `analysis/` scripts
Make sure there is a `calc/` folder in the main directory, as that's where the output data from the simulation is put. Afterwards, the plot.py script can be used to visualize the results.

Some nice commands:

    ./run.sh animate # short-hand for `make; ./2dsynch animate; py3 analysis/plot.py animate`
    ./run.sh lossmap
    ./run.sh startdist

### 2. Input for SixTrack
A special tool has been designed to output data directly to SixTrack.

Compile:

    make export
    
Use by:

    ./export coll . 1 10000 linexp 1000

The output is found at `./1.txt`

Arguments explained:
1. Exports for SixTrack collimation version or non-collimation version: `coll` or `nocoll`
2. Directory file path for output (relative or absolute). Write `.` for the current directory. Output files are named `1.txt`, `2.txt`, ... for as many files as is requested in argument 3. Normal rules for `relative/paths` and `/absolute/paths` apply. 
3. Number of output files
4. Number of particles per file (should be an even number)
5. Longitudinal distribution. Can be any of the following: `linexp`, `lin`, `const`
    - `linexp`: Linear distribution inside separatrix, exponential outside. This is the final estimated distribution that was found.
    - `lin`: Linear distribution from the separatrix to the collimator.
    - `const`: Constant distribution from the separatrix to the collimator.
6. Number of turns to simulate


# Analysis tools
Quick overview of the different files:

**Note: The term 'lossmap' here is a little outdated and wrongly used. It is the terminology I used for simulated BLM signals.**

- `plot.py` : Main tool for visualising output from the Toy model. Can make movies/distribution plots/loss time profiles. Some useful commands:
    - `py3 analysis/plot.py dist <path/to/distribution>`
    - `py3 analysis/plot.py animate`
    - `py3 analysis/plot.py ham-dist`
    - `py3 analysis/plot.py e-dist`
- `batch_plot.py` : Same role as `plot.py` but for LXPLUS batches.
- `phasespace.py` : Data structure for managing phase space data for all particles. It reads `calc/startdist.dat`, `calc/enddist.dat`, `calc/particles.dat`, `calc/lines.dat` depending on what is being studied.
- `lossmap.py` : Data structure for keeping track of what turn a particle hits a collimator (the data structure has been renamed to `CHitMap`). Reads `calc/TCP_IR3.ch` and `calc/TCPc_IR7.ch`.
- `lhccomp.py` : Compares Toy model with LHC aggregate. Aggregate fill generated from `../lhcstat/` needs to exist to work.
- `ramp.py` : Reads the energy ramp from file.
- `sixtrack_oml.py` : Reads in SixTrack batches, and does plots similar to `plot.py` (e.g., plotting distributions, time profiles). Does not do a proper SixTrack loss map.
- `settings.py` and `settings.json` : Contains some project settings. 
