## ----setup, include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, echo = TRUE, collapse = TRUE, fig.align = "center")
knitr::purl("main.Rmd", documentation = 1)


## ----gui, echo=FALSE, fig.cap="**QGIS**, an example of Graphical User Interface (GUI) software", out.width="80%"----
knitr::include_graphics("images/QGIS.png")


## ----cli, echo=FALSE, fig.cap="**R**, an example of Command Line Interface (CLI) software", out.width="80%"----
knitr::include_graphics("images/R.png")


## ---- echo=FALSE, fig.cap="Books on Spatial Data Analysis with R", out.width="100%"----
knitr::include_graphics("images/books.svg")


## ----buffer, echo=FALSE, results="hide", message=FALSE, fig.cap="Buffer function", fig.width=6, fig.height=2.5, out.width="100%", warning=FALSE----
library(sf)
nafot = st_read("data/nafot.shp")
p = nafot[nafot$name_eng == "Be'er Sheva", ]
p = st_geometry(p)
opar = par(mfrow=c(1,4), mar = c(0, 0, 1, 0))

plot(p %>% st_buffer(100000), main = "10 km buffer", border = NA)
plot(p, add = TRUE)
plot(p %>% st_buffer(10000) %>% st_difference(p), col = "lightgrey", add = TRUE)

plot(p %>% st_buffer(100000), main = "25 km buffer", border = NA)
plot(p, add = TRUE)
plot(p %>% st_buffer(25000) %>% st_difference(p), col = "lightgrey", add = TRUE)

plot(p %>% st_buffer(100000), main = "50 km buffer", border = NA)
plot(p, add = TRUE)
plot(p %>% st_buffer(50000) %>% st_difference(p), col = "lightgrey", add = TRUE)

plot(p %>% st_buffer(100000), main = "100 km buffer")
plot(p, add = TRUE)
plot(p %>% st_buffer(100000) %>% st_difference(p), col = "lightgrey", add = TRUE)

par(opar)


## ----gstat, echo=FALSE, results="hide", message=FALSE, fig.cap="Predicted Zinc concentration, using Ordinary Kriging", out.width="70%"----
library(gstat)
library(automap)
library(stars)

# Prepare data
data(meuse)
data(meuse.riv)
coordinates(meuse) = ~ x + y
data(meuse.grid)
gridded(meuse.grid) = ~ x + y
grid = st_as_stars(meuse.grid)

# Predict
f = log(zinc) ~ 1
v = autofitVariogram(f, meuse)
g = gstat(formula = log(zinc) ~ 1, model = v$var_model, data = meuse)
predicted = predict(g, grid)

# River color
col = col2rgb("lightblue")
col = col / 255
col = rgb(col[1], col[2], col[3], 0.5)

# Plot
plot(predicted, col = rev(hcl.colors(11, "Spectral")), reset = FALSE, key.pos = 4, main = NA)
polygon(meuse.riv, asp = 1, col = col)
plot(meuse, pch = 1, cex = log(meuse$zinc) / 5, add = TRUE)


## ----spdep, echo=FALSE, results="hide", message=FALSE, warning=FALSE, fig.cap="Neighbors list based on regions with contiguous boundaries", fig.width=7, fig.height=3.5, out.width="100%"----
# From '?poly2nb'
library(spdep)
nc = st_read(system.file("shape/nc.shp", package = "sf"))
nc = as(nc, "Spatial")
nc$rate = nc$SID79 / nc$BIR79
nb = poly2nb(nc)
opar = par(mar = rep(0, 4))
plot(nc, border = "grey")
plot(nb, coordinates(nc), add = TRUE, col = "black")
par(opar)


## ----spatstat, echo=FALSE, results="hide", message=FALSE, warning=FALSE, fig.cap="Distance map for the Biological Cells point pattern dataset", fig.width=4, fig.height=4, out.width="50%"----
library(spatstat)
data(cells)
U = distmap(cells)
opar = par(mar = rep(0, 4))
contour(U, main = "")
plot(cells, add = TRUE, col = "red", pch = 3)
par(opar)


## ---- eval=FALSE---------------------------------------------------------------------
## # library(sf)
## # library(RPostgreSQL)
## 
## # con = dbConnect(
## #   PostgreSQL(),
## #   dbname = "gisdb",
## #   host = "159.89.13.241",
## #   port = 5432,
## #   user = "geobgu",
## #   password = "*******"
## # )


## ---- eval=FALSE---------------------------------------------------------------------
## # q = "SELECT name_lat, geometry FROM plants LIMIT 3;"
## # st_read(con, query = q)


## ---- echo=FALSE---------------------------------------------------------------------
# library(sf)
# library(RPostgreSQL)
# con = dbConnect(
#   PostgreSQL(),
#   dbname = "gisdb",
#   host = "159.89.13.241",
#   port = 5432,
#   user = "geobgu",
#   password = "geobgu1"
# )
# q = "SELECT name_lat, geometry FROM plants LIMIT 3;"
# st_read(con, query = q)


## ---- include=FALSE------------------------------------------------------------------
# dbDisconnect(con)


## ----sf-r-journal, echo=FALSE, fig.cap="Pebesma, 2018, The R Journal (https://journal.r-project.org/archive/2018-1/)", out.width="100%"----
knitr::include_graphics("images/lesson_07_paper1.png")


## ----geometry-types, echo=FALSE, fig.cap="Simple Feature types (see also: https://r-spatial.github.io/sf/articles/sf1.html)", fig.width=6.8, fig.height=4, out.width="80%", warning=FALSE, message=FALSE----
library(sf)
library(units)
point = st_as_sfc("POINT (30 10)")[[1]]
linestring = st_as_sfc("LINESTRING (30 10, 10 30, 40 40)")[[1]]
polygon = st_as_sfc("POLYGON ((35 10, 45 45, 15 40, 10 20, 35 10),(20 30, 35 35, 30 20, 20 30))")[[1]]
multipoint = st_as_sfc("MULTIPOINT ((10 40), (40 30), (20 20), (30 10))")[[1]]
multilinestring = st_as_sfc("MULTILINESTRING ((10 10, 20 20, 10 40),(40 40, 30 30, 40 20, 30 10))")[[1]]
multipolygon = st_as_sfc("MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)),((20 35, 10 30, 10 10, 30 5, 45 20, 20 35),(30 20, 20 15, 20 25, 30 20)))")[[1]]
geometrycollection = st_as_sfc("GEOMETRYCOLLECTION (POLYGON((30 20, 45 40, 10 40, 30 20)),LINESTRING (10 10, 20 20, 10 30),POINT (40 20))")[[1]]
pol = st_as_sfc("POLYGON((30 20, 45 40, 10 40, 30 20))")[[1]]
l = st_as_sfc("LINESTRING (10 10, 20 20, 10 30)")[[1]]
p = st_as_sfc("POINT (40 20)")[[1]]
opar = par(mfrow = c(2, 4), mar = c(1,1,1,1))
plot(point, main = "POINT", col = "blue", cex = 1.8, lwd = 2)
plot(linestring, main = "LINESTRING", col = "blue", lwd = 2)
plot(polygon, main = "POLYGON", border = "blue", col = "#0000FF33", lwd = 2)
plot(1, type="n", axes=F, xlab="", ylab="")
plot(multipoint, main = "MULTIPOINT", col = "blue", cex = 1.8, lwd = 2)
plot(multilinestring, main = "MULTILINESTRING", col = "blue", lwd = 2)
plot(multipolygon, main = "MULTIPOLYGON", border = "blue", col = "#0000FF33", lwd = 2)
plot(geometrycollection, main = "GEOMETRYCOLLECTION", col = NA, border = NA, lwd = 2)
plot(pol, border = "blue", col = "#0000FF33", add = TRUE, lwd = 2)
plot(l, col = "blue", add = TRUE, lwd = 2)
plot(p, col = "blue", add = TRUE, cex = 1.8, lwd = 2)
par(opar)


## ----sf-dependencies, echo=FALSE, fig.cap="`sf` package dependencies (https://github.com/edzer/rstudio_conf)", out.width="100%"----
knitr::include_graphics("images/sf_deps.png")


## ----nc-geometry-column, echo=FALSE, fig.cap="Structure of an `sf` object (https://cran.r-project.org/web/packages/sf/vignettes/sf1.html)", out.width="100%"----
knitr::include_graphics("images/sf.png")


## ----nc-plot, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Visualization of the `sf` object shown in Figure \\@ref(fig:nc-geometry-column)", out.width="100%", fig.width=8, fig.height=3----
library(sf)
nc = st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc = nc[1:3, 9:15]

library(ggplot2)
nc$id = 1:nrow(nc)
nc1 = reshape2::melt(st_drop_geometry(nc), id.vars = "id")
nc1 = merge(nc1, nc, all.x = TRUE)
nc1 = st_sf(nc1)
ctr = st_centroid(nc1)
ctr$x = st_coordinates(ctr)[, 1]
ctr$y = st_coordinates(ctr)[, 2]

opar = par(mar = rep(0, 4))
ggplot() +
  geom_sf(data = nc1) +
  geom_text(data = ctr, aes(x = x, y = y, label = value)) +
  facet_wrap(~ variable) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = 'transparent')
  )
par(opar)


## ---- include=FALSE------------------------------------------------------------------
lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))


## ------------------------------------------------------------------------------------
library(sf)


## ------------------------------------------------------------------------------------
nafot = st_read("data/nafot.shp")


## ------------------------------------------------------------------------------------
nafot


## ------------------------------------------------------------------------------------
class(nafot)


## ------------------------------------------------------------------------------------
st_geometry(nafot)


## ------------------------------------------------------------------------------------
st_drop_geometry(nafot)


## ----simple-plot, fig.cap="The `nafot` layer", fig.width=4, fig.height=6, out.width="50%"----
plot(nafot, key.width = lcm(4))


## ---- eval=FALSE---------------------------------------------------------------------
## plot(st_geometry(nafot))


## ----simple-plot-geometry, echo=FALSE, fig.cap="The `nafot` layer", fig.width=4, fig.height=6, out.width="50%"----
opar = par(mar = rep(0, 4))
plot(st_geometry(nafot))
par(opar)


## ----geographic-vs-projected, echo=FALSE, fig.cap="US counties in WGS84 and US Atlas projections", out.width="100%"----
knitr::include_graphics("images/geographic-vs-projected-2.png")


## ------------------------------------------------------------------------------------
st_crs(nafot)


## ------------------------------------------------------------------------------------
nafot_wgs84 = st_transform(nafot, 4326)
nafot_wgs84


## ---- fig.cap="Nafot in UTM and WGS84 coordinate reference systems", fig.height=7, fig.width=4.5, out.width="40%", fig.show="hold"----
plot(st_geometry(nafot), main = "UTM", axes = TRUE)
plot(st_geometry(nafot_wgs84), main = "WGS84", axes = TRUE)


## ------------------------------------------------------------------------------------
rail = st_read("data/RAIL_STRATEGIC.shp")


## ------------------------------------------------------------------------------------
rail


## ------------------------------------------------------------------------------------
rail = st_transform(rail, st_crs(nafot))


## ----plot-rail, fig.cap="The `rail` layer", fig.width=6, fig.height=4, out.width="100%", warning=FALSE----
plot(rail)


## ---- eval=FALSE---------------------------------------------------------------------
## plot(st_geometry(nafot), border = "grey")
## plot(st_geometry(rail), add = TRUE)


## ----nafot-rail, echo=FALSE, fig.cap="The `nafot` and `rail` geometries", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
opar = par(mar = rep(0, 4))
plot(st_geometry(nafot), border = "grey")
plot(st_geometry(rail), add = TRUE)
par(opar)


## ------------------------------------------------------------------------------------
rail = rail[rail$ISACTIVE == "פעיל", ]
rail


## ------------------------------------------------------------------------------------
rail$segment_id = 1:nrow(rail)
rail = rail["segment_id"]
rail


## ------------------------------------------------------------------------------------
nafot1 = nafot[rail, ]
nafot1


## ---- eval=FALSE---------------------------------------------------------------------
## plot(st_geometry(nafot), border = "grey50")
## plot(st_geometry(nafot1), border = "grey50", col = "grey90", add = TRUE)
## plot(st_geometry(rail), add = TRUE)


## ----nafot-subset, echo=FALSE, fig.cap="The `nafot` and `rail` geometries", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
opar = par(mar = rep(0, 4))
plot(st_geometry(nafot), border = "grey50")
plot(st_geometry(nafot1), border = "grey50", col = "grey90", add = TRUE)
plot(st_geometry(rail), add = TRUE)
par(opar)


## ------------------------------------------------------------------------------------
nafot$area = st_area(nafot)
nafot$area[1:3]


## ------------------------------------------------------------------------------------
class(nafot$area)


## ------------------------------------------------------------------------------------
library(units)
nafot$area = set_units(nafot$area, "km^2")
nafot$area[1:3]


## ---- fig.cap="Calculated `area` attribute", fig.width=4, fig.height=6, out.width="50%"----
plot(nafot[, "area"])


## ------------------------------------------------------------------------------------
int = st_intersects(nafot1, nafot1, sparse = FALSE)


## ------------------------------------------------------------------------------------
int


## ---- fig.cap="Intersection relations between `nafot` features", out.width="80%"-----
int1 = apply(int, 2, rev)
int1 = t(int1)
image(int1, col = c("lightgrey", "red"), asp = 1, axes = FALSE)
axis(3, at = seq(0, 1, 1/(nrow(int1)-1)), labels = nafot1$name_eng, las = 2, lwd = 0, lwd.ticks = 1, cex.axis = 0.75)
axis(2, at = seq(0, 1, 1/(nrow(int1)-1)), labels = rev(nafot1$name_eng), las = 1, pos = -0.046, lwd = 0, lwd.ticks = 1, cex.axis = 0.75)


## ---- echo=FALSE, fig.cap="Geometry-generating operations on individual layers", fig.height=4, fig.width=6, out.width="100%"----
set.seed(1)
x = st_multipoint(matrix(runif(10), ncol = 2))
x = st_buffer(st_sfc(lapply(1:3, function(x) st_point(c(x,x)))), 0.2 * 1:3)
opar = par(mfrow=c(2,3), mar = rep(1,4))
plot(x, border = '#ff333388')
plot(st_centroid(x), add = TRUE, pch = 3)
title("st_centroid")
plot(x, border = '#ff333388')
plot(st_buffer(x, dist = 0.1), add = TRUE, pch = 3)
plot(st_buffer(x, dist = 0.2), add = TRUE, pch = 3)
plot(st_buffer(x, dist = 0.3), add = TRUE, pch = 3)
plot(st_buffer(x, dist = 0.4), add = TRUE, pch = 3)
plot(st_buffer(x, dist = 0.5), add = TRUE, pch = 3)
title("st_buffer")
s = split(x, 1:3)
s = lapply(s, st_sample, size = 5)
s = lapply(s, st_combine)
s = do.call(c, s)
plot(x, border = '#ff333388')
plot(s, add = TRUE, pch = 3)
title("st_sample")
plot(s, col = '#ff333388', pch = 3)
plot(st_convex_hull(s), add = TRUE, pch = 3)
title("st_convex_hull")
s = st_union(s)
v = st_voronoi(s)
plot(s, col = '#ff333388', pch = 3)
plot(v, col = NA, border = 1, axes = FALSE, add = TRUE)
title("st_voronoi")
par(opar)


## ---- warning=FALSE------------------------------------------------------------------
nafot_ctr = st_centroid(nafot)


## ---- eval=FALSE---------------------------------------------------------------------
## plot(st_geometry(nafot), border = "grey")
## plot(st_geometry(nafot_ctr), col = "red", pch = 3, add = TRUE)


## ----nafot-ctr, echo=FALSE, fig.cap="State centroids", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
opar = par(mar = rep(0, 4))
plot(st_geometry(nafot), border = "grey")
plot(st_geometry(nafot_ctr), col = "red", pch = 3, add = TRUE)
par(opar)


## ----geometry-generating-pairs, echo=FALSE, fig.cap="Geometry-generating operations on pairs of layers", fig.height = 3.5, fig.width = 6, out.width="80%"----
x = st_point(c(0, 0))
x = st_buffer(x, 0.5)
y = st_point(c(0.5, 0))
y = st_buffer(y, 0.5)

xy = c(x, y)

opar = par(mfrow=c(2,3), mar = rep(1,4))

plot(xy, border = NA)
plot(x, add = TRUE, col = "#ff333388")
plot(y, add = TRUE, col = "#33ff3388")
title("x: red, y: green")

plot(xy, border = "grey")
plot(st_intersection(x, y), col = "lightblue", add = TRUE)
title("intersection(x, y)")

plot(xy, border = "grey")
plot(st_difference(x, y), col = "lightblue", add = TRUE)
title("difference(x, y)")

plot(xy, border = NA)

plot(xy, border = "grey")
plot(st_sym_difference(x, y), col = "lightblue", add = TRUE)
title("sym_difference(x, y)")

plot(xy, border = "grey")
plot(st_union(x, y), col = "lightblue", add = TRUE)
title("union(x, y)")

par(opar)


## ---- message=FALSE, warning=FALSE---------------------------------------------------
rail_int = st_intersection(rail, nafot)
rail_int


## ---- fig.cap="Intersection result", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
plot(rail_int[, "name_eng"], lwd = 3, key.width = lcm(4), reset = FALSE)
plot(st_geometry(nafot), border = "lightgrey", add = TRUE)


## ------------------------------------------------------------------------------------
class(rail_int$geometry)


## ------------------------------------------------------------------------------------
rail_int = st_cast(rail_int, "MULTILINESTRING")
rail_int


## ------------------------------------------------------------------------------------
rail_int$length = st_length(rail_int)
rail_int$length = set_units(rail_int$length, km)
rail_int


## ---- message=FALSE------------------------------------------------------------------
rail_int = st_drop_geometry(rail_int) 
rail_int = aggregate(rail_int["length"], rail_int["name_eng"], sum)


## ------------------------------------------------------------------------------------
head(rail_int)


## ------------------------------------------------------------------------------------
nafot = merge(nafot, rail_int, by = "name_eng", all.x = TRUE)


## ------------------------------------------------------------------------------------
nafot


## ------------------------------------------------------------------------------------
nafot$length[is.na(nafot$length)] = 0


## ---- fig.cap="Total railway length per Nafa", out.width="50%", warning=FALSE, fig.width=4, fig.height=6----
plot(nafot[, "length"])


## ------------------------------------------------------------------------------------
nafot$density = nafot$length / nafot$area


## ------------------------------------------------------------------------------------
nafot = nafot[order(nafot$density, decreasing = TRUE), ]
nafot


## ----density-1, fig.cap="railway density per state", out.width="50%", warning=FALSE, fig.width=4, fig.height=6----
plot(nafot[, "density"])


## ----density-2, fig.cap="Nafot layer with calculated attributes", warning=FALSE------
plot(nafot)


## ------------------------------------------------------------------------------------
stat = st_read("data/statisticalareas_demography2018.gdb")


## ------------------------------------------------------------------------------------
vars = colnames(stat)
vars


## ------------------------------------------------------------------------------------
vars = vars[grepl("age_", vars)]
vars


## ------------------------------------------------------------------------------------
stat = stat[vars]


## ------------------------------------------------------------------------------------
stat = st_transform(stat, st_crs(nafot))


## ----stat-layer, fig.width=9, fig.height=5, fig.cap="Population estimates in the `stat` layer", out.width="100%"----
plot(stat, max.plot = 16)


## ----stat-one-attr, fig.width=4, fig.height=6, fig.cap="The `age_10_14` attribute in `stat`, and the `nafot` layer", out.width="50%"----
plot(stat["age_10_14"], pal = hcl.colors(12, "Reds", rev = TRUE), border = "black", lwd = 0.07, reset = FALSE)
plot(st_geometry(nafot), border = "grey", add = TRUE)


## ------------------------------------------------------------------------------------
stat[is.na(stat)] = 0


## ----stat-one-attr2, fig.width=4, fig.height=6, fig.cap="The `age_10_14` attribute in `stat`, and the `nafot` layer", out.width="50%"----
plot(stat["age_10_14"], pal = hcl.colors(12, "Reds", rev = TRUE), border = "black", lwd = 0.07, reset = FALSE)
plot(st_geometry(nafot), border = "grey", add = TRUE)


## ---- eval=FALSE---------------------------------------------------------------------
## st_interpolate_aw(stat, nafot, extensive = TRUE)


## ---- error=TRUE---------------------------------------------------------------------
x = st_interpolate_aw(stat, nafot, extensive = TRUE)


## ------------------------------------------------------------------------------------
as.data.frame(table(st_geometry_type(stat)))


## ------------------------------------------------------------------------------------
stat = st_cast(stat, "MULTIPOLYGON")
as.data.frame(table(st_geometry_type(stat)))


## ---- error=TRUE---------------------------------------------------------------------
x = st_interpolate_aw(stat, nafot, extensive = TRUE)


## ------------------------------------------------------------------------------------
stat = st_make_valid(stat)


## ------------------------------------------------------------------------------------
x = st_interpolate_aw(stat, nafot, extensive = TRUE)
x$Group.1 = NULL


## ---- fig.width=9, fig.height=5, out.width="100%"------------------------------------
plot(x, max.plot = 16)


## ------------------------------------------------------------------------------------
dat = st_drop_geometry(x)
rownames(dat) = nafot$name_eng
dat[, 1:6]


## ------------------------------------------------------------------------------------
dat = sweep(dat, 1, rowSums(dat), "/")
dat[, 1:6]


## ---- fig.width=7, fig.height=4, out.width="100%"------------------------------------
d = dist(dat)
hc = hclust(d, "average")
k = 3
groups = cutree(hc, k = k)
plot(hc)
rect.hclust(hc, k = k)


## ------------------------------------------------------------------------------------
library(reshape2)

# Reshape
dat2 = dat
dat2$group = groups
dat2$name = rownames(dat2)
dat2 = melt(dat2, id.vars = c("group", "name"))
dat2$variable = gsub("age_", "", dat2$variable)

# Set age groups order
x = strsplit(dat2$variable, "_")
x = sapply(x, "[", 1)
x = as.numeric(x)
x = dat2$variable[order(x)]
x = unique(x)
dat2$variable = factor(dat2$variable, levels = x)

# Set cluster labels
x = split(dat2$name, dat2$group)
x = sapply(x, "[", 1)
dat2$group = factor(dat2$group, levels = 1:3, labels = x)

head(dat2)


## ----ggplot-1, fig.cap="Nafa clusters (lines)", fig.width=8, fig.height=4, out.width="100%"----
library(ggplot2)
ggplot(dat2, aes(x = variable, y = value, group = name)) +
    geom_line() +
    scale_x_discrete("Age group") +
    scale_y_continuous("Proportion") +
    facet_wrap(~ group) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )


## ----ggplot-2, fig.cap="Nafa clusters (columns)", fig.width=8, fig.height=4, out.width="100%"----
ggplot(dat2, aes(x = name, y = value, fill = variable, group = name)) +
    geom_col(colour = "black") +
    scale_x_discrete("Nafa") +
    scale_y_continuous("Proportion") +
    scale_fill_manual("Age group", values = hcl.colors(length(levels(dat2$variable)), "Dark2")) +
    facet_grid(. ~ group, scales = "free", space = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )

