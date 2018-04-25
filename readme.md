# HiMAP: High-resolution Microbial Analysis Pipeline

## Installation from source

### Linux

For Linux, we need to make sure you have development version of OpenSSL installed, so e.g. for Ubuntu run from terminal:
```
sudo apt-get -y install libssl-dev libcurl4-openssl-dev
```
Then to install R, we need to add the repository, authentication key and install it:
```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
sudo apt-get update
sudo apt-get -y install r-base
```

In R, install devtools package, then use `install_github` to install the package from source:

```R
install.packages('devtools')
devtools::install_github('benjjneb/dada2')
# Select 1 when prompted to install biocLite
devtools::install_github('taolonglab/himap', auth_token='99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6')
```

This may take up to about 30 minutes to install on a fresh R installation.

### macOS

First, install Xcode command line tools. To do this, open a macOS Terminal, then type
```
xcode-select --install
```

```R
options('install.packages.compile.from.source'='never')
install.packages('devtools')

# DADA2 pre-requisites
source('http://bioconductor.org/biocLite.R')
biocLite('GenomeInfoDbData', suppressUpdates=T)

# DADA2
devtools::install_github('benjjneb/dada2')
# Select 1 when prompted to install biocLite

# Finally install HiMAP
devtools::install_github('taolonglab/himap', auth_token='99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6')
```

If you don't have the Xcode command line tools, you will be prompted to install them. After that's completed, re-run `devtools::install_github` commands.


## Tutorials

* [Human gut microbiome data tutorial](tutorial.ipynb)
* [mock community tutorial](tutorial_mock.ipynb)

## Notes

By default, multithreading is enabled and supported only on Linux and macOS, through R package `parallel`.


