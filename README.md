# clas12_analysis

Quick test to show dst2root is working and output a few histograms to a root file.

Prerequisites:
* [cern root](https://root.cern/install/)

## cpp
A c++ example can be found in the cpp folder.

To build:

```bash
git clone https://github.com/tylern4/clas12_analysis.git
mkdir -p clas12_analysis/build
cd clas12_analysis/build
cmake ..
make
```

To run:
```bash
BEAM_E=7.5 NUM_THREADS=4 ./clas12_analysis output.root /path/to/input/*.root
```

```bash
BEAM_E=10.6 NUM_THREADS=4 ./clas12_yields yields.csv /path/to/input/*.root
```

There are three folders under `src`, `lib` has all the `cpp` files which will have most of the class definitions, logic, and functions. The `include` folder holds all the header files with the definitions in it. The executables are found in `exe` and are usually the same, get input from the command line and use a class (defined in `include`/`lib`) to start running a specific process. 

## python

This uses a small wrapper library I wrote around [boost histogram](https://boost-histogram.readthedocs.io/en/latest/) to plot the way I wanted in python. [nicks_plot_utils](https://github.com/tylern4/nicks_plot_utils) can be installed with `pip install nicks-plot-utils.`

To run:

```
BEAM_E=10.6 python main.py '/path/to/data/*.root'
```

Not a full list but some needed packages from pip/conda.

```
boost-histogram
nicks-plot-utils
pandas
scipy
matplotlib
numpy
lmfit
```

And then you would need root installed with the python libraries installed to use `TLorentzVectors`.