# HiMAP: High-resolution Microbial Analysis Pipeline

## Installation

For Linux, make sure you have development version of OpenSSL installed, so e.g. for Ubuntu run from terminal:
```
sudo apt install libssl-dev
```
Add these lines to the end of `/etc/apt/sources.list`:
```
deb http://cloud.r-project.org/bin/linux/ubuntu xenial/
deb http://cloud.r-project.org/ trusty-backports main restricted universe
```
Install R:
```
sudo apt install r-base
```
(ignore any warnings about packages that cannot be authenticated).


In R, install devtools package, then use `install_github` to install the package from source:

```R
install.packages('devtools')
devtools::install_github('benjjneb/dada2')
# Select 1 when prompted to install biocLite
devtools::install_github('taolonglab/himap', auth_token='99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6')
```

This may take up to about 30 minutes to install on a fresh R installation.


## Tutorials

* [Human gut microbiome data tutorial](tutorial.ipynb)
* [mock community tutorial](tutorial_mock.ipynb)

## Notes

By default, multithreading is enabled and supported only on Linux and macOS, through R
