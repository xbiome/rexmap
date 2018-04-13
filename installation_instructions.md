# HiMAP installation from source

Now install HiMAP from github (need OAuth token for private repo):
```R
install.packages('devtools')
devtools::install_github('benjjneb/dada2')
devtools::install_github('taolonglab/himap', auth_token='99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6')
```
This should install dependent packages automatically.
