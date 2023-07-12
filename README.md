# DLOPy_v1.0

Python code for determining orientation of OBS horizontals written by Adrian. K. Doran and Gabi Laske (see [Doran & Laske, 2017](https://github.com/jbrussell/DLOPy_v1.0/blob/master/README/BSSA-2016165.1.pdf)). It automatically searches for events and downloads data from the IRIS-DMC to compute the orientations.

By default, it assumes H2 is 90 degrees CW from H1 with Z pointing up (left handed system). It uses global dispersion maps to predict the Rayleigh-wave arrival window. The cross-correlations are preformed in seven frequency bands ranging from 10 to 40 mHz (25-100 s).

This version is allows for automated looping over many stations.

**run_orientations.py** - Calculate station orientations. Specify network, station, event parameters at the top of the file. Data are automatically downloaded from IRIS.

**run_orientations_local.py** - Same, but this version reads local miniseed data.

