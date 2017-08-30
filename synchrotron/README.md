Small project aimed to help understand the longitudinal motion of particles inside a synchrotron accelerator. 

# Compiling

The code has been tested on macOS Sierra and on LXPLUS. To compile it the first time, simply write `make` in the main directory.

Options:
- It is possible to run the code multi-threaded. Go into the `makefile` and change `MULTI_THREADING=1`. Note that [Threading Building Blocks (TBB)](https://www.threadingbuildingblocks.org) is needed for multi-threading. Installing it on LXPLUS is a little bit tricky, so in general I recommend running the single threaded version there. Installing TBB on macOS is easy using [MacPorts](https://www.macports.org), and it should be possible to use `apt-get` on Linux.

# Use the code
There are 2 main use cases for the tool:

### Simulation meant to be analysed using `analysis/` scripts.
Make sure there is a calc/ folder in the main directory, as that's where the output data from the simulation is put. Afterwards, the plot.py script can be used to visualize the results.

Some nice commands:

    ./run.sh animate # short-hand for `make; ./2dsynch animate; py3 analysis/plot.py animate`
    ./run.sh lossmap
    ./run.sh startdist

### Input for SixTrack
A special tool has been designed to output data directly to SixTrack.

Compile:

    make export
    
Use by:

    ./export coll . 1 10000 linexp 1000

The output is found at `./1.txt`.

Arguments explained:
1. Exports for SixTrack collimation version or non-collimation version: `coll` or `nocoll`
2. Directory file path for output (relative or absolute). Write `.` for the current directory. Output files are named `1.txt`, `2.txt`, ... for as many files as is requested in argument 3.
    E.g.: 'relative/path', '/absolute/path'. 
3. Number of output files
4. Number of particles per file (should be an even number)
5. Longitudinal distribution. Can be any of the following: `linexp`, `lin`, `const`
    - `linexp`: Linear distribution inside separatrix, exponential outside. This is the final estimated distribution that was found.
    - `lin`: Linear distribution from the separatrix to the collimator.
    - `const`: Constant distribution from the separatrix to the collimator.
6. Number of turns to simulate
