## ---- eval=FALSE-------------------------------------------------------
## install.packages("sf")


## ---- message=TRUE, warning=TRUE---------------------------------------
library(sf)


## ---- eval=FALSE-------------------------------------------------------
## setwd("C:\\Data\\CSB")


## ----------------------------------------------------------------------
url = "https://github.com/michaeldorman/R-Spatial-Workshop-at-CBS-2021/raw/main/data.zip"
download.file(url, "data.zip")
unzip("data.zip")


## ----------------------------------------------------------------------
nafot = st_read("data/nafot.shp")


## ----------------------------------------------------------------------
rail = st_read("data/RAIL_STRATEGIC.shp", options = "ENCODING=UTF-8")


## ----------------------------------------------------------------------
stat = st_read("data/statisticalareas_demography2018.gdb", options = "ENCODING=UTF-8")


## ---- eval=FALSE-------------------------------------------------------
## nafot = st_read("nafot.shp")
## rail = st_read("RAIL_STRATEGIC.shp")
## stat = st_read("statisticalareas_demography2018.gdb")


## ----nafot, fig.cap='The `nafot` layer'--------------------------------
plot(nafot, key.width = lcm(4))


## ----rail, fig.cap='The `rail` layer', fig.width=6, fig.height=4, out.width="100%"----
plot(rail)


## ----stat, fig.cap='The `stat` layer'----------------------------------
plot(stat)


## ----------------------------------------------------------------------
rail = rail[rail$ISACTIVE == "פעיל", ]


## ----rail2, fig.cap='The `rail` layer, subset of active railway lines', fig.width=7, fig.height=3, out.width="100%"----
plot(rail)


## ----encoding, echo=FALSE, fig.cap="Reopening an R script with the UTF-8 encoding in RStudio", out.width="100%"----
knitr::include_graphics("images/encoding.png")

