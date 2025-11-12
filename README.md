# Pathopath

Utilities to build patient pathways and compute contact networks from ward-level movements,
with quick plotting and a simple Shiny app.

## Install from source (local folder)
```r
install.packages(".", repos = NULL, type = "source", dependencies = TRUE)
library(pathopath)
```

## Try the included datasets
```r
chrf <- read.csv(system.file("extdata","CHRF_2021.csv", package = "pathopath"))
demo <- read.csv(system.file("extdata","patients_demo.csv", package = "pathopath"))
```

## Shiny app
```R
source(system.file("scripts","Pathopath_shiny.R", package = "pathopath"))
```
