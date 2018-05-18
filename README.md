## Using K-means Analysis to Classify Synthetic Spectral Profiles from Numerical Simulations of the Solar Atmosphere

## Introduction
We classify the synthetic profiles of MG II k & h in our dataset (544005 profiles) using k-means analysis. We show the locations of the k-means for certain cluster numbers as well as the spectral profile for k-mean label corresponding to those location.This code allows to read RH (in netCDF format) profiles. It allows to selects specific lines, interpolate in an uniform spectra and perform k-mean analysis. 

## Installation
Works on MacOS, Linux, and Windows. 

    [] Download the Python 3.6 version (can also download 2.7 version if necessary)
    [] Clone the latest version of helita from Github: git clone https://github.com/jamiehuang00/K-Means-IRIS
    [] Change helita/helita/io/__init__.py to read:
        __all__ = ["crispex", "fio", "lp"] #, "ncdf", "sdf"]
        from . import crispex
        from . import lp
        from . import ncdf

Requires Python 3.0 or higher.

## Usage


## Contributing
If you would like to improve the code, your help is welcomed. To make a contribution, go to the right of the Branch menu, and click New pull request.

## Credits
Juan Martinez-Sykora, Alberto Sainz-Dalda
