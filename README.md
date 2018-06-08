## Using K-means Analysis to Classify Synthetic Spectral Profiles from Numerical Simulations of the Solar Atmosphere

## Introduction
We classify the synthetic profiles of MG II k & h in our dataset (544005 profiles) using k-means analysis. We show the locations of the k-means for certain cluster numbers as well as the spectral profile for k-mean label corresponding to those location.This code allows to read RH (in netCDF format) profiles. It allows to selects specific lines, interpolate in an uniform spectra and perform k-mean analysis. 

## Installation
Works on MacOS, Linux, and Windows.

1) Download the Python 3.6 version (can also download 2.7 version if necessary) 

2) Clone the latest version of helita from Github: git clone https://github.com/jamiehuang00/K-Means-IRIS and download to desktop

3) Use Terminal to compile and run the code

Requires Python 3.0 or higher.

## Usage
In order to run the code, enter this on terminal or Jupyter Notebook: 

python
***

ipython

>>> import numpy as np

>>> import matplotlib.pyplot as plt

>>> from sklearn.cluster import KMeans

>>> from sklearn.cluster import MiniBatchKMeans

>>> import pickle

>>> import imp

>>> import helita

>>> cd ~/helita/helita/sim

>>> from helita import kmeans as km

>>> __init__

>>> interp()

>>> time_import(tm, inertia, t0)

>>> fit(t0)

>>> create_km_map()

>>> create_k_means_maps(self, wvl)

>>> create_spectral_map(self, i3, wvl, ax)

## Contributing
If you would like to improve the code or report a bug, your help is welcomed. 
Here are the steps:

1) Fork the repository.
2) Develop and test code changes.
3) Verify that tests pass successfully.
4) Start discussion or give feedback 
5) Look over your changes in the diffs on the Compare page, make sure they’re what you want to submit.
6) Push to your fork repository
7) Go to the right of the Branch menu
8) Select the master branch, and click New pull request.

## Credits
Juan Martinez-Sykora, Alberto Sainz-Dalda
