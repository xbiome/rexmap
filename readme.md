# HiMAP: High-resolution Microbial Analysis Pipeline

## Installation

### Linux

For Linux, make sure you have development version of OpenSSL installed, so e.g. for Ubuntu run from terminal:
```
sudo apt-get -y install libssl-dev libcurl4-openssl-dev
```
Add these lines to the end of `/etc/apt/sources.list`, e.g. `sudo nano /etc/apt/sources.list`:
```
deb http://cloud.r-project.org/bin/linux/ubuntu xenial/
```
Install the latest version of R:
```
sudo apt-get --allow-unauthenticated update
sudo apt-get -y --allow-unauthenticated install r-base
```
(ignore any warnings about packages that cannot be authenticated). At this point you can install RStudio.


In R, install devtools package, then use `install_github` to install the package from source:

```R
install.packages('devtools')
devtools::install_github('benjjneb/dada2')
# Select 1 when prompted to install biocLite
devtools::install_github('taolonglab/himap', auth_token='99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6')
```

This may take up to about 30 minutes to install on a fresh R installation.

### macOS



## Tutorials

* [Human gut microbiome data tutorial](tutorial.ipynb)
* [mock community tutorial](tutorial_mock.ipynb)

## Notes

By default, multithreading is enabled and supported only on Linux and macOS, through R
