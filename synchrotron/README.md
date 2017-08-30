Small project aimed to help understand the longitudinal motion of particles inside a synchrotron accelerator. 

# Compiling

The code has been tested on macOS Sierra and on LXPLUS. To compile it the first time, simply write `make` in the main directory.

Options:
- It is possible to run the code multi-threaded. Go into the `makefile` and change `MULTI_THREADING=1`. Note that [Threading Building Blocks (TBB)](https://www.threadingbuildingblocks.org) is needed for multi-threading. Installing it on LXPLUS is a little bit tricky, so in general I recommend running the single threaded version there. Installing TBB on macOS is easy using [MacPorts](https://www.macports.org), and it should be possible to use `apt-get` on Linux.

# Use the Toy model
There are 2 main use cases for the Toy model:

### 1. Simulation meant to be analysed using `analysis/` scripts
To get started, do the following steps:

1. Create a `calc/` directory at the install location (same directory as this file is in)
2. (Optional) If you want your simulations cached, also create a directory called `simulations/cache/`
3. Open `settings.json`:
    - Edit `DATA_DIR` to the absolute path where the `calc/` directory is located
    - Edit `RESOURCE_DIR` to the absolute path where the `resource/` directory is located 

Now you should be good to go! Some nice commands:

    ./run.sh animate # short-hand for `make; ./2dsynch animate; py3 analysis/plot.py animate`
    ./run.sh lossmap
    ./run.sh startdist

**Note: The term 'lossmap' here is a little outdated and wrongly used. It is the terminology I used for simulated BLM signals.**

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
Dependencies:
- Python 3.5
- matplotlib
- numpy
- scipy


Quick overview of the main files:
- `plot.py` : Main tool for visualising output from the Toy model. Can make movies/distribution plots/loss time profiles. 
    Some useful commands:

      python3.5 analysis/plot dist <path/to/distribution>   # Plots the phase space of the given distribution
      python3.5 analysis/plot animate                       # Animates the phase space of the distribution found in `calc/particles.dat`
      python3.5 analysis/plot ham-dist                      # Plots a histogram over the `calc/startdist.dat` in Hamiltonian units
      python3.5 analysis/plot e-dist                        # Same as above, but for ΔE units instead

- `batch_plot.py` : Same role as `plot.py` but for LXPLUS batches.
- `phasespace.py` : Data structure for managing phase space data for all particles. It reads `calc/startdist.dat`, `calc/enddist.dat`, `calc/particles.dat`, `calc/lines.dat` depending on what is being studied.
- `lossmap.py` : Data structure for keeping track of what turn a particle hits a collimator (the data structure has been renamed to `CHitMap`). Reads `calc/TCP_IR3.ch` and `calc/TCPc_IR7.ch`.
- `lhccomp.py` : Compares Toy model with LHC aggregate. Aggregate fill generated from `../lhcstat/` needs to exist to work.
