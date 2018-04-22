# HiMAP: High-resolution Microbial Analysis Pipeline

## Installation

In R, install devtools package, then use `install_github` to install the package from source:

```R
install.packages('devtools')
devtools::install_github('taolonglab/himap', auth_token='99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6')
```

Every dependency should install automatically. If DADA2 fails install the latest version from source, before HiMAP:
```R
devtools::install_github('benjjneb/dada2')
```

## Tutorials

* [Human gut microbiome data tutorial](tutorial.ipynb)
* [mock community tutorial](tutorial_mock.ipynb)

## Notes

By default, multithreading is enabled and supported only on Linux and macOS, through R