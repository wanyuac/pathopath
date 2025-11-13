# Pathopath

Utilities to build patient pathways and compute contact networks from ward-level movements,
with quick plotting and a simple Shiny app.

## Installation
```r
install.packages(".", repos = NULL, type = "source", dependencies = TRUE)
library(pathopath)
```

## Prepare input files

Please ensure input files do not have any missing values ("" or NA).

* `movements.tsv` from user's records of locations and time. The column `Location` stands for location accessions, which are independent to any specific location level (Hospital, Building, Floor, Ward, Unit, Bed, _etc_).
  * `Time_in` and `Time_out`: time stamps following [ISO 8601](https://www.iso.org/iso-8601-date-and-time-format.html). Pathopath current only supports dates (YYYY-MM-DD).

* `genotypes.tsv`: currently not supported.

## Try the included datasets

```r
chrf <- read.csv(system.file("extdata","CHRF_2021.csv", package = "pathopath"))
demo <- read.csv(system.file("extdata","patients_demo.csv", package = "pathopath"))
```

## Shiny app
```R
source(system.file("scripts","Pathopath_shiny.R", package = "pathopath"))
```

## Strengths of pathopath

* Versability: Support multiple location levels (Hospital, Building, Floor, Unit, Ward, Room, Bed, *etc*) that can be specified by users.
* Generality: Incorporation of patients and inanimate subjects.
