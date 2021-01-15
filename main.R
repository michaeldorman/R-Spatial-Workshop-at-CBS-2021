## ------------------------------------------------------
library(sf)


## ------------------------------------------------------
nafot = st_read("data/nafot.shp")


## ------------------------------------------------------
nafot


## ------------------------------------------------------
class(nafot)


## ------------------------------------------------------
st_geometry(nafot)


## ------------------------------------------------------
st_drop_geometry(nafot)


## ----simple-plot, fig.cap="The `nafot` layer", fig.width=4, fig.height=6, out.width="50%"----
plot(nafot, key.width = lcm(4))


## ----simple-plot-geometry, echo=FALSE, fig.cap="The `nafot` layer", fig.width=4, fig.height=6, out.width="50%"----
plot(st_geometry(nafot))


## ------------------------------------------------------
st_crs(nafot)


## ------------------------------------------------------
nafot_wgs84 = st_transform(nafot, 4326)


## ------------------------------------------------------
nafot_wgs84


## ----nafot-wgs84, fig.cap="Nafot in UTM and WGS84 coordinate reference systems", fig.height=7, fig.width=4.5, out.width="40%", fig.show="hold"----
plot(st_geometry(nafot), main = "UTM", axes = TRUE)
plot(st_geometry(nafot_wgs84), main = "WGS84", axes = TRUE)


## ------------------------------------------------------
rail = st_read("data/RAIL_STRATEGIC.shp", options = "ENCODING=UTF-8")


## ------------------------------------------------------
rail


## ------------------------------------------------------
rail = st_transform(rail, st_crs(nafot))


## ----plot-rail, fig.cap="The `rail` layer", fig.width=6, fig.height=4, out.width="100%", warning=FALSE----
plot(rail)


## ----plot-rail-one-attr, fig.cap="One attribute of the `rail` layer", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
plot(rail["ISACTIVE"], key.width = lcm(3))


## ----nafot-rail, echo=FALSE, fig.cap="The `nafot` and `rail` geometries", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
plot(st_geometry(nafot), border = "grey")
plot(st_geometry(rail), col = "red", add = TRUE)


## ------------------------------------------------------
rail = rail[rail$ISACTIVE == "פעיל", ]
rail


## ------------------------------------------------------
rail$segment_id = 1:nrow(rail)
rail = rail["segment_id"]
rail


## ----rail-active, fig.cap="Subset of active railway lines", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
plot(rail)


## ------------------------------------------------------
nafot1 = nafot[rail, ]
nafot1


## ----nafot-subset, echo=FALSE, fig.cap="The `nafot` and `rail` geometries", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
plot(st_geometry(nafot), border = "grey")
plot(st_geometry(nafot1), border = "grey", col = "grey90", add = TRUE)
plot(st_geometry(rail), col = "red", add = TRUE)


## ------------------------------------------------------
nafot$area = st_area(nafot)


## ------------------------------------------------------
nafot


## ------------------------------------------------------
class(nafot$area)


## ------------------------------------------------------
library(units)


## ------------------------------------------------------
nafot$area = set_units(nafot$area, "km^2")


## ------------------------------------------------------
nafot


## ----nafot-area, fig.cap="Calculated `area` attribute", fig.width=4, fig.height=6, out.width="50%"----
plot(nafot["area"])


## ------------------------------------------------------
int = st_intersects(nafot1, nafot1, sparse = FALSE)


## ------------------------------------------------------
int


## ----st-intersects-matrix, fig.cap="Intersection relations between `nafot` features", out.width="80%"----
int1 = apply(int, 2, rev)
int1 = t(int1)
image(int1, col = c("lightgrey", "red"), asp = 1, axes = FALSE)
axis(3, at = seq(0, 1, 1/(nrow(int1)-1)), labels = nafot1$name_eng, las = 2, lwd = 0, lwd.ticks = 1, cex.axis = 0.75)
axis(2, at = seq(0, 1, 1/(nrow(int1)-1)), labels = rev(nafot1$name_eng), las = 1, pos = -0.046, lwd = 0, lwd.ticks = 1, cex.axis = 0.75)


## ---- warning=FALSE------------------------------------
nafot_ctr = st_centroid(nafot)


## ----nafot-ctr, echo=FALSE, fig.cap='"Nafot" centroids', fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
plot(st_geometry(nafot), border = "grey")
plot(st_geometry(nafot_ctr), col = "red", pch = 3, add = TRUE)


## ---- message=FALSE, warning=FALSE---------------------
rail = st_union(rail)
rail


## ------------------------------------------------------
rail = st_sf(geometry = rail)
rail


## ---- message=FALSE, warning=FALSE---------------------
rail = st_intersection(rail, nafot)


## ------------------------------------------------------
rail


## ----nafa-rail-intersection, fig.cap="Intersection result", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
plot(rail["name_eng"], lwd = 3, key.width = lcm(4), reset = FALSE)
plot(st_geometry(nafot), border = "lightgrey", add = TRUE)


## ------------------------------------------------------
rail$length = st_length(rail)
rail


## ------------------------------------------------------
rail$length = set_units(rail$length, "km")
rail


## ----nafa-rail-lengths, fig.cap="Intersection result", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----
plot(rail["name_eng"], lwd = 3, key.width = lcm(4), reset = FALSE)
plot(st_geometry(nafot), border = "lightgrey", add = TRUE)
text(st_coordinates(st_centroid(rail)), as.character(round(rail$length)))


## ------------------------------------------------------
rail = st_drop_geometry(rail)
rail


## ------------------------------------------------------
rail = rail[c("name_eng", "length")]
rail


## ------------------------------------------------------
nafot = merge(nafot, rail, by = "name_eng", all.x = TRUE)


## ------------------------------------------------------
nafot


## ------------------------------------------------------
nafot$length[is.na(nafot$length)] = 0
nafot


## ----length-per-nafa, fig.cap="Total railway length per Nafa", out.width="50%", warning=FALSE, fig.width=4, fig.height=6----
plot(nafot[, "length"])


## ------------------------------------------------------
nafot$density = nafot$length / nafot$area


## ------------------------------------------------------
nafot = nafot[order(nafot$density, decreasing = TRUE), ]
nafot


## ----density-1, fig.cap='Railway density per "Nafa"', out.width="50%", warning=FALSE, fig.width=4, fig.height=6----
plot(nafot["density"])


## ----density-2, fig.cap="Nafot layer with calculated attributes", warning=FALSE----
plot(nafot)


## ------------------------------------------------------
stat = st_read("data/statisticalareas_demography2018.gdb", options = "ENCODING=UTF-8")


## ------------------------------------------------------
vars = colnames(stat)
vars


## ------------------------------------------------------
vars = vars[grepl("age_", vars)]
vars
stat = stat[vars]


## ----stat-layer, fig.width=9, fig.height=5, fig.cap="Population estimates in the `stat` layer", out.width="100%"----
# plot(stat, max.plot = 16)


## ------------------------------------------------------
stat[is.na(stat)] = 0


## ----stat-layer-filled, fig.width=9, fig.height=5, fig.cap="Population estimates in the `stat` layer", out.width="100%"----
# plot(stat, max.plot = 16)


## ------------------------------------------------------
stat = st_transform(stat, st_crs(nafot))


## ----stat-one-attr, fig.width=4, fig.height=6, fig.cap="The `age_10_14` attribute in `stat`, and the `nafot` layer", out.width="50%"----
plot(stat["age_10_14"], pal = hcl.colors(12, "Reds", rev = TRUE), border = "black", lwd = 0.07, reset = FALSE)
plot(st_geometry(nafot), lwd = 0.5, add = TRUE)


## ---- warning=FALSE, error=TRUE------------------------
# x = st_interpolate_aw(stat, nafot, extensive = TRUE)


## ------------------------------------------------------
as.data.frame(table(st_geometry_type(stat)))


## ------------------------------------------------------
stat = st_cast(stat, "MULTIPOLYGON")


## ------------------------------------------------------
as.data.frame(table(st_geometry_type(stat)))


## ------------------------------------------------------
stat = st_make_valid(stat)


## ---- warning=FALSE------------------------------------
x = st_interpolate_aw(stat, nafot, extensive = TRUE)


## ----aw-nafa, fig.cap='Area-weighted counts per age group, by "Nafa"', fig.width=9, fig.height=5, out.width="100%"----
plot(x, max.plot = 16)


## ------------------------------------------------------
dat = st_drop_geometry(x)
dat[, 1:5]


## ------------------------------------------------------
rownames(dat) = nafot$name_eng
dat[, 1:5]


## ------------------------------------------------------
dat = sweep(dat, 1, rowSums(dat), "/")
dat[, 1:5]


## ------------------------------------------------------
d = dist(dat)


## ------------------------------------------------------
hc = hclust(d, "average")


## ------------------------------------------------------
k = 3
groups = cutree(hc, k = k)
groups


## ----hc-dendrogram, fig.cap="Hierarchical dendrogram of age structure per Nafa", fig.width=7, fig.height=5, out.width="100%"----
plot(hc)
rect.hclust(hc, k = k)


## ----hc-map, fig.cap="Age structure clusters on a map", fig.width=4, fig.height=6, out.width="50%"----
nafot$group = factor(groups)
plot(nafot["group"])


## ------------------------------------------------------
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


## ----ggplot-1, fig.cap="Distribution of age structures in each cluster of `nafot`", fig.width=12, fig.height=4, out.width="100%"----
library(ggplot2)
library(ggrepel)
ggplot(dat2, aes(x = variable, y = value, group = name, colour = name)) +
    geom_line() +
    geom_point() +
    geom_label_repel(data = dat2[dat2$variable == "75_up", ], aes(label = name)) +
    scale_x_discrete("Age group") +
    scale_y_continuous("Proportion") +
    scale_colour_discrete(guide = FALSE) +
    facet_wrap(~ group) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
    )

