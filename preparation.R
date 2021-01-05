## ----setup, include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, echo = TRUE, collapse = TRUE, fig.align = "center")
knitr::purl("preparation.Rmd", documentation = 1)


## ---- eval=FALSE---------------------------------------------------------------------
## install.packages("sf")


## ---- message=TRUE-------------------------------------------------------------------
library(sf)


## ---- eval=FALSE---------------------------------------------------------------------
## setwd("C:\\Data\\CSB")


## ------------------------------------------------------------------------------------
url = "https://github.com/michaeldorman/R-Spatial-Workshop-at-CBS-2021/raw/main/data.zip"
download.file(url, "data.zip")
unzip("data.zip")


## ------------------------------------------------------------------------------------
nafot = st_read("data/nafot.shp")
rail = st_read("data/RAIL_STRATEGIC.shp", options = "ENCODING=UTF-8")
stat = st_read("data/statisticalareas_demography2018.gdb", options = "ENCODING=UTF-8")


## ---- eval=FALSE---------------------------------------------------------------------
## nafot = st_read("nafot.shp")
## rail = st_read("RAIL_STRATEGIC.shp")
## stat = st_read("statisticalareas_demography2018.gdb")


## ------------------------------------------------------------------------------------
plot(nafot, key.width = lcm(4))


## ------------------------------------------------------------------------------------
plot(rail)


## ------------------------------------------------------------------------------------
plot(stat)


## ------------------------------------------------------------------------------------
rail = rail[rail$ISACTIVE == "פעיל", ]


## ------------------------------------------------------------------------------------
plot(rail)


## ----encoding, echo=FALSE, fig.cap="Reopening an R script with the UTF-8 encoding in RStudio", out.width="100%"----
knitr::include_graphics("images/encoding.png")
