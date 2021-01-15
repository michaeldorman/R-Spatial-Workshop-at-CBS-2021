## ---- eval=FALSE---------------------------------------
## install.packages("sf")


## ---- message=TRUE, warning=TRUE-----------------------
library(sf)

## ---- eval=FALSE---------------------------------------
## url = "https://github.com/michaeldorman/R-Spatial-Workshop-at-CBS-2021/raw/main/data.zip"
## download.file(url, "data.zip")
## unzip("data.zip")


## ------------------------------------------------------
nafot = st_read("data/nafot.shp")


## ------------------------------------------------------
rail = st_read("data/RAIL_STRATEGIC.shp", options = "ENCODING=UTF-8")


## ------------------------------------------------------
stat = st_read("data/statisticalareas_demography2018.gdb", options = "ENCODING=UTF-8")

## ----nafot, fig.cap='The `nafot` layer'----------------
plot(nafot, key.width = lcm(4))


## ----rail, fig.cap='The `rail` layer', fig.width=6, fig.height=4, out.width="100%"----
plot(rail)


## ----stat, fig.cap='The `stat` layer'------------------
plot(stat)


## ------------------------------------------------------
rail = rail[rail$ISACTIVE == "פעיל", ]


## ----rail2, fig.cap='The `rail` layer, subset of active railway lines', fig.width=7, fig.height=3, out.width="100%"----
plot(rail)

