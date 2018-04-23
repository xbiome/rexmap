# HiMAP: High-resolution Microbial Analysis Pipeline

## Installation

For Ubuntu, make sure you have development version of OpenSSL installed, so run from terminal:
```
sudo apt install libssl-dev
```

In R, install devtools package, then use `install_github` to install the package from source:

```R
install.packages('devtools')
devtools::install_github('benjjneb/dada2')
# Select 1 when prompted to install biocLite
devtools::install_github('taolonglab/himap', auth_token='99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6')
```

This may take about 15 minutes to install.


## Tutorials

* [Human gut microbiome data tutorial](tutorial.ipynb)
* [mock community tutorial](tutorial_mock.ipynb)

## Notes

By default, multithreading is enabled and supported only on Linux and macOS, through R
