#' ## Reading vector data
#' 
#' We move on to practical examples. First, let's see how creating and working with `sf` layers works in practice. 
#' 
#' Function `st_read` can be used to read vector layers into `sf` data structures. Before doing anything else we need to load the `sf` package:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(sf)

#' 
#' Then, we can use `st_read` to import a vector into the R environment. In this case, we are reading the `nafot.shp` Shapefile:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot = st_read("data/nafot.shp")

#' 
#' Printing the object gives a summary of its properties, and the values of the (first 10) features. In this case, we see that the `nafot` layer has `r nrow(nafot)` features, and one attribute called `name_eng`:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot

#' 
#' The `class` function, applied on an `sf` layer, returns both `"sf"` and `"data.frame"` class names:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class(nafot)

#' 
#' As mentioned above, a layer (geometry+attributes) is represented by an `sf` object. The second printed value, `data.frame`, implies that `sf` is, in fact, an extension of the `data.frame` class, inheriting many of its methods. Indeed, we will see that many functions in R work exactly the same way on `sf` and `data.frame` objects.
#' 
#' If we want just the *geometric* part, it can be extracted with `st_geometry`, resulting in an object of class `sfc`, i.e., a geometry column (Table \@ref(tab:structures-package-sf)):
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
st_geometry(nafot)

#' 
#' Conversely, If we want just the *non-geometric* part (the "attributes"), it can be extracted with `st_drop_geometry` which returns a `data.frame`:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
st_drop_geometry(nafot)

#' 
#' The `plot` function is a quick way to see the spatial arrangment and attribute values in an `sf` layer. For example, the `nafot` layer attribute can be plotted as follows (Figure \@ref(fig:simple-plot)):
#' 
## ----simple-plot, fig.cap="The `nafot` layer", fig.width=4, fig.height=6, out.width="50%"----------------------------------------------------------------------------------------------------------------
plot(nafot, key.width = lcm(4))

#' 
#' Alternatively, we can plot the geometry only, without attributes, as follows (Figure \@ref(fig:simple-plot-geometry)):
#' 
## ---- eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## plot(st_geometry(nafot))

#' 
## ----simple-plot-geometry, echo=FALSE, fig.cap="The `nafot` layer", fig.width=4, fig.height=6, out.width="50%"-------------------------------------------------------------------------------------------
plot(st_geometry(nafot))

#' 
#' ## Coordinate Reference Systems (CRS)
#' 
#' ### What are CRS?
#' 
#' A **Coordinate Reference System (CRS)** defines how the coordinates in our geometries relate to the surface of the Earth. There are two main types of CRS:
#' 
#' * **Geographic**---longitude and latitude, in degrees
#' * **Projected**---implying flat surface, usually in units of true distance (e.g., meters)
#' 
#' For example, Figure \@ref(fig:geographic-vs-projected) shows the same polygonal layer (U.S. counties) in two different projection. On the left, the county layer is in the WGS84 geographic projection. Indeed, we can see that the axes are given in degrees of longitude and latitude. For example, we can see how the nothern border of U.S. follows the 49° latitude line. On the right, the same layer is displayed in the US National Atlas projection, where units are arbitrary but reflect true distance (meters). For example, the distance between every two consecutive grid lines is 1,000,000 meters or 1,000 kilometers.
#' 
## ----geographic-vs-projected, echo=FALSE, fig.cap="US counties in WGS84 and US Atlas projections", out.width="100%"--------------------------------------------------------------------------------------

#' 
#' The CRS of a given `sf` layer can be obtained using function `st_crs`. The CRS information is returned in the WKT format:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
st_crs(nafot)

#' 
#' ### Vector layer reprojection {#vector-layer-reprojection}
#' 
#' **Reprojection** is the transformation of geometry coordinates, from one CRS to another. It is an important part of spatial analysis workflow, since we often need to:
#'     
#' * Transform several layers into the same projection, so that they can be displayed one on top of the other (e.g., Figure \@ref(fig:nafot-rail)) or so that they can be subject to a spatial operator (e.g., Figure \@ref(fig:geometry-generating-pairs))
#' * Switch between geographic and projected CRS
#' 
#' A vector layer can be reprojected with `st_transform`. The `st_transform` function has two important parameters: 
#' 
#' * `x`---The layer to be reprojected
#' * `crs`---The *target* CRS
#' 
#' The `crs` can be specified in one of four ways:
#' 
#' * An **EPSG** code (e.g., `4326`)
#' * A **PROJ4** string (e.g., `"+proj=longlat +datum=WGS84 +no_defs"`)
#' * A **WKT** string
#' * A `crs` object of another layer, as returned by `st_crs`
#' 
#' For example, the following expression reprojects the `nafot` layer to the geographic WGS84 CRS, using its EPSG code (`4326`):
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot_wgs84 = st_transform(nafot, 4326)

#' 
#' Printing the layer summary demostrates that the coordinates and CRS definition have changed:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot_wgs84

#' 
#' The two layers, `nafot` and `nafot_wgs84`, are shown in Figure \@ref(fig:nafot-wgs84):
#' 
## ----nafot-wgs84, fig.cap="Nafot in UTM and WGS84 coordinate reference systems", fig.height=7, fig.width=4.5, out.width="40%", fig.show="hold"-----------------------------------------------------------
plot(st_geometry(nafot), main = "UTM", axes = TRUE)
plot(st_geometry(nafot_wgs84), main = "WGS84", axes = TRUE)

#' 
#' # Geoprocessing functions
#' 
#' ## Reading layers into R
#' 
#' For the next series of examples, we read a second layer called `RAIL_STRATEGIC.shp`. The `RAIL_STRATEGIC.shp` layer contains railway lines in Israel. Using `st_read`, we import the layer and create an `sf` object called `rail`:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail = st_read("data/RAIL_STRATEGIC.shp", options = "ENCODING=UTF-8")

#' 
#' (Note that the `ENCODING` option is only needed when operating system default encoding is not the same as the layer encoding.)
#' 
#' Here is a summary of the `rail` layer:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail

#' 
#' ## Reprojection
#' 
#' For any type of spatial analysis, we usually need all input layers to be in the same CRS. For this purpose, we will reproject the `rail` layer to the CRS of the `nafot` layer:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail = st_transform(rail, st_crs(nafot))

#' 
#' ## Basic plotting
#' 
#' We can plot the `rail` layer on its own, as shown in Figure \@ref(fig:plot-rail), to examine its attributes:
#' 
## ----plot-rail, fig.cap="The `rail` layer", fig.width=6, fig.height=4, out.width="100%", warning=FALSE---------------------------------------------------------------------------------------------------
plot(rail)

#' 
#' We can also plot the `nafot` and `rail` geometries together, to examine their arrangement (Figure \@ref(fig:nafot-rail)). The second expression uses `add=TRUE` to add the geometries on top of the existing plot, instead of initializing a new one:
#' 
## ---- eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## plot(st_geometry(nafot), border = "grey")
## plot(st_geometry(rail), col = "red", add = TRUE)

#' 
## ----nafot-rail, echo=FALSE, fig.cap="The `nafot` and `rail` geometries", fig.width=4, fig.height=6, out.width="50%", warning=FALSE----------------------------------------------------------------------
plot(st_geometry(nafot), border = "grey")
plot(st_geometry(rail), col = "red", add = TRUE)

#' 
#' ## Subsetting
#' 
#' ### Non-spatial
#' 
#' Subsetting (filtering) of features in an `sf` vector layer is done in exactly the same way as filtering rows in a `data.frame`. For example, the following expression filters the `rail` layer to keep only those railway lines which are active: 
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail = rail[rail$ISACTIVE == "פעיל", ]
rail

#' 
#' We can also subset the columns (attributes) we need. For example, in the following expressions we create an ID attribute called `segment_id`, and remove all other attributes:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail$segment_id = 1:nrow(rail)
rail = rail["segment_id"]
rail

#' 
#' The modified layer is shown in Figure \@ref(fig:rail-active):
#' 
## ----rail-active, fig.cap="Subset of active railway lines", fig.width=4, fig.height=6, out.width="50%", warning=FALSE------------------------------------------------------------------------------------
plot(rail)

#' 
#' ### Spatial
#' 
#' We can also subset features according to whether they *intersect* with another layer, using the latter as an index. For example, the following expression creates a subset of `nafot`, named `nafot1`, with only those features intersecting the `rail` layer:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot1 = nafot[rail, ]
nafot1

#' 
#' Figure \@ref(fig:nafot-subset) shows the `nafot1` subset (in grey fill), the complete `nafot` layer (in grey outline), and the railway lines layer (in red):
#' 
## ---- eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## plot(st_geometry(nafot), border = "grey")
## plot(st_geometry(nafot1), border = "grey", col = "grey90", add = TRUE)
## plot(st_geometry(rail), col = "red", add = TRUE)

#' 
## ----nafot-subset, echo=FALSE, fig.cap="The `nafot` and `rail` geometries", fig.width=4, fig.height=6, out.width="50%", warning=FALSE--------------------------------------------------------------------
plot(st_geometry(nafot), border = "grey")
plot(st_geometry(nafot1), border = "grey", col = "grey90", add = TRUE)
plot(st_geometry(rail), col = "red", add = TRUE)

#' 
#' ## Geometric calculations
#' 
#' Geometric operations on vector layers can conceptually be divided into three groups according to their output:
#' 
#' * **Numeric** values (Section \@ref(numeric-geometric-calculations))---Functions that summarize geometrical properties of:
#'     * A *single* layer---e.g., area, length
#'     * A *pair* of layers---e.g., distance
#' * **Logical** values (Section \@ref(logical-geometric-calculations))---Functions that evaluate whether a certain condition holds true, regarding:
#'     * A *single* layer---e.g., geometry is valid
#'     * A *pair* of layers---e.g., feature A intersects feature B
#' * **Spatial** layers (Section \@ref(spatial-geometric-calculations))---Functions that create a new layer based on:
#'     * A *single* layer---e.g., centroid, buffer
#'     * A *pair* of layers---e.g., intersection area
#' 
#' ### Numeric {#numeric-geometric-calculations}
#' 
#' There are several functions to calculate *numeric* geometric properties of vector layers in package `sf`:
#' 
#' * `st_length`
#' * `st_area`
#' * `st_distance`
#' * `st_bbox`
#' * `st_dimension`
#' 
#' For example, we can calculate the area of each feature in the `nafot` layer (i.e. each state) using `st_area`:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot$area = st_area(nafot)

#' 
#' The `nafot` layer now has a new attribute with the polygon areas:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot

#' 
#' The resulting column is an object of class `units`. This is basically a numeric vector associated with units of measurement:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class(nafot$area)

#' 
#' We can convert measurements to different units with `set_units`, from package `units`:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(units)
nafot$area = set_units(nafot$area, "km^2")

#' 
#' The modified `area` column can be observed in the layer summary:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot

#' 
#' Figure \@ref(fig:nafot-area) shows a map of the calculated area sizes:
#' 
## ----nafot-area, fig.cap="Calculated `area` attribute", fig.width=4, fig.height=6, out.width="50%"-------------------------------------------------------------------------------------------------------
plot(nafot["area"])

#' 
#' ### Logical {#logical-geometric-calculations}
#' 
#' Given two layers, `x` and `y`, the following *logical* geometric functions check whether each feature in `x` maintains the specified relation with each feature in `y`:
#' 
#' * `st_intersects`
#' * `st_disjoint`
#' * `st_touches`
#' * `st_crosses`
#' * `st_within`
#' * `st_contains`
#' * `st_overlaps`
#' * `st_covers`
#' * `st_covered_by`
#' * `st_equals`
#' * `st_equals_exact`
#' 
#' When specifying `sparse=FALSE` the functions return a logical `matrix`. Each element `i,j` in the `matrix` is `TRUE` when `f(x[i], y[j])` is `TRUE`. For example, this creates a matrix of *intersection* relations between nafot:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int = st_intersects(nafot1, nafot1, sparse = FALSE)

#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int

#' 
#' The following code section visualizes the `matrix` (Figure \@ref(fig:st-intersects-matrix)):
#' 
## ----st-intersects-matrix, fig.cap="Intersection relations between `nafot` features", out.width="80%"----------------------------------------------------------------------------------------------------
int1 = apply(int, 2, rev)
int1 = t(int1)
image(int1, col = c("lightgrey", "red"), asp = 1, axes = FALSE)
axis(3, at = seq(0, 1, 1/(nrow(int1)-1)), labels = nafot1$name_eng, las = 2, lwd = 0, lwd.ticks = 1, cex.axis = 0.75)
axis(2, at = seq(0, 1, 1/(nrow(int1)-1)), labels = rev(nafot1$name_eng), las = 1, pos = -0.046, lwd = 0, lwd.ticks = 1, cex.axis = 0.75)

#' 
#' ### Spatial {#spatial-geometric-calculations}
#' 
#' `sf` provides common *geometry-generating* functions applicable to individual geometries, such as:
#' 
#' * `st_centroid`
#' * `st_buffer`
#' * `st_sample`
#' * `st_convex_hull`
#' * `st_voronoi`
#' 
## ---- echo=FALSE, fig.cap="Geometry-generating operations on individual layers", fig.height=4, fig.width=6, out.width="80%"------------------------------------------------------------------------------

#' 
#' For example, the following expression uses `st_centroid` to create a layer of "Nafa" centroids:
#' 
## ---- warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot_ctr = st_centroid(nafot)

#' 
#' They can be plotted as follows, the result is shown in Figure \@ref(fig:nafot-ctr):
#' 
## ---- eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## plot(st_geometry(nafot), border = "grey")
## plot(st_geometry(nafot_ctr), col = "red", pch = 3, add = TRUE)

#' 
## ----nafot-ctr, echo=FALSE, fig.cap="State centroids", fig.width=4, fig.height=6, out.width="50%", warning=FALSE-----------------------------------------------------------------------------------------
plot(st_geometry(nafot), border = "grey")
plot(st_geometry(nafot_ctr), col = "red", pch = 3, add = TRUE)

#' 
#' Other geometry-generating functions work on *pairs* of input geometries (Figure \@ref(fig:geometry-generating-pairs)):
#' 
#' * `st_intersection`
#' * `st_difference`
#' * `st_sym_difference`
#' * `st_union`
#' 
## ----geometry-generating-pairs, echo=FALSE, fig.cap="Geometry-generating operations on pairs of layers", fig.height = 3.5, fig.width = 6, out.width="80%"------------------------------------------------
#' 

#' Now we move on to examples of two specific tasks:
#' 
#' * Calculating line density per polygon (Section \@ref(example-1))
#' * Calculating heirarchical clusters (Section \@ref(example-2))
#' 
#' # Example #1: Rail density {#example-1}
#' 
#' ## Splitting lines by polygons
#' 
#' In this example, our goal is to calculate total rail density (length per $km^2$) per "Nafa" (Figure \@ref(fig:density-2)). First, we will use `st_intersection` to 'split' the `rail` layer by "Nafa":
#' 
## ---- message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail_int = st_intersection(rail, nafot)

#' 
#' The result is a new line layer of `rail` segments, split by "Nafa" borders and including the `name_eng` attribute of the corresponding "Nafa":
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail_int

#' 
#' The layer is visualized in Figure \@ref(fig:nafa-rail-intersection):
#' 
## ----nafa-rail-intersection, fig.cap="Intersection result", fig.width=4, fig.height=6, out.width="50%", warning=FALSE------------------------------------------------------------------------------------
plot(rail_int["name_eng"], lwd = 3, key.width = lcm(4), reset = FALSE)
plot(st_geometry(nafot), border = "lightgrey", add = TRUE)

#' 
#' ## Line length
#' 
#' Next, we will calculate line **length** of each rail segment, using function `st_length`. We the result into a railway length attribute called `length`:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail_int$length = st_length(rail_int)
rail_int

#' 
#' For convenience, we will also convert the values from $m$ to $km$:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail_int$length = set_units(rail_int$length, "km")
rail_int

#' 
#' ## Aggregation
#' 
#' Since we are interested in total length per "Nafa", we will aggregate the attribute table of `rail_int` by `state`, to find the sum of `length` values. 
#' 
#' First, we discard the geometries, keeping just the attributes:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail_int = st_drop_geometry(rail_int) 
head(rail_int)

#' 
#' Then, we use `aggregate` to calculate total segments length per "Nafa":
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail_int = aggregate(rail_int["length"], rail_int["name_eng"], sum)

#' 
#' The result is a `data.frame` with total length of railway tracks, per "Nafa":
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rail_int

#' 
#' ## Join layer with table
#' 
#' To attach the track lengths to the `nafot` layer, we use a non-spatial *join*:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot = merge(nafot, rail_int, by = "name_eng", all.x = TRUE)

#' 
#' As a result of the join operation, the `nafot` layer now has a new `length` attribute:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot

#' 
#' Length of `NA` implies there are no railways in that polygon. These `NA` values should be replaced with zero:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot$length[is.na(nafot$length)] = 0
nafot

#' 
#' The calculated total lengths per "Nafa" are shown in Figure \@ref(fig:length-per-nafa):
#' 
## ----length-per-nafa, fig.cap="Total railway length per Nafa", out.width="50%", warning=FALSE, fig.width=4, fig.height=6---------------------------------------------------------------------------------
plot(nafot[, "length"])

#' 
#' ## Calculating density
#' 
#' Now that we know the total railway length and the area, per "Nafa", we can divide total railway length by area. This gives us railway **density** per "Nafa":
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot$density = nafot$length / nafot$area

#' 
#' Sorting the features lets us clearly see which "Nafa" has the highest densities of railway lines:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nafot = nafot[order(nafot$density, decreasing = TRUE), ]
nafot

#' 
#' Plotting the layer shows the new `density` attribute (Figures \@ref(fig:density-1) and \@ref(fig:density-2)):
#' 
## ----density-1, fig.cap="railway density per state", out.width="50%", warning=FALSE, fig.width=4, fig.height=6-------------------------------------------------------------------------------------------
plot(nafot["density"])

#' 
## ----density-2, fig.cap="Nafot layer with calculated attributes", warning=FALSE--------------------------------------------------------------------------------------------------------------------------
plot(nafot)

#' 
#' # Example #2: Population patterns {#example-2}
#' 
#' ## Reading statistical areas
#' 
#' In the second example, we are going to examine the dissimilarity between "Nafot" in terms of their age structure. The demographic data come at the statistical area level, which we are going to aggregate to the "Nafa" level. 
#' 
#' First, we read the statistical areas Shapefile:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
stat = st_read("data/statisticalareas_demography2018.gdb", options = "ENCODING=UTF-8")

#' 
#' ## Subsetting
#' 
#' The layer has numerous columns:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
vars = colnames(stat)
vars

#' 
#' but, in this case, we are only interested in the population estimates per age group:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
vars = vars[grepl("age_", vars)]
vars

#' 
#' We will retain only the latter the attributes:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
stat = stat[vars]

#' 
#' ## Reprojecting
#' 
#' Again, we need to make sure both layers are in the same CRS:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
stat = st_transform(stat, st_crs(nafot))

#' 
#' The resulting subset of the `stat` layer is shown in Figure \@ref(fig:stat-layer).
#' 
## ----stat-layer, fig.width=9, fig.height=5, fig.cap="Population estimates in the `stat` layer", out.width="100%"-----------------------------------------------------------------------------------------
plot(stat, max.plot = 16)

#' 
#' Figure \@ref(fig:stat-one-attr) shows one of the attributes in `stat` with the `nafot` layer on top:
#' 
## ----stat-one-attr, fig.width=4, fig.height=6, fig.cap="The `age_10_14` attribute in `stat`, and the `nafot` layer", out.width="50%"---------------------------------------------------------------------
plot(stat["age_10_14"], pal = hcl.colors(12, "Reds", rev = TRUE), border = "black", lwd = 0.07, reset = FALSE)
plot(st_geometry(nafot), border = "grey", add = TRUE)

#' 
#' ## Fill missing data
#' 
#' One thing we may notice in Figure \@ref(fig:stat-one-attr), is that many of the statistical areas have `NA` values, meaning zero (rather than unknown) population size. For all practical purposes these should be replaced with "true" zero:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
stat[is.na(stat)] = 0

#' 
#' The modification is demonstrated in Figure \@ref(fig:stat-one-attr2). Now we are ready to "transfer" the demographic estimates from the statistical area level, to the "Nafa" level.
#' 
## ----stat-one-attr2, fig.width=4, fig.height=6, fig.cap="The `age_10_14` attribute in `stat`, and the `nafot` layer", out.width="50%"--------------------------------------------------------------------
plot(stat["age_10_14"], pal = hcl.colors(12, "Reds", rev = TRUE), border = "black", lwd = 0.07, reset = FALSE)
plot(st_geometry(nafot), border = "grey", add = TRUE)

#' 
#' ## Area weighted sum
#' 
#' Now, in theory, we are ready to calculate the weighted sum of population per "Nafa", using function `st_interpolate_aw`. The 
#' function requires three parameters:
#' 
#' * `x`---Input layer (with values to transfer)
#' * `to`---Target geometries (where to transfer values)
#' * `extensive`---Spatially extensive, like population (`TRUE`), or spatially intensive, like population density (`FALSE`)
#' 
## ---- eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## st_interpolate_aw(stat, nafot, extensive = TRUE)

#' 
#' ## Standardizing geometries
#' 
#' Now, in theory, we are ready to calculate the weighted sum of population per "Nafa", using function `st_interpolate_aw`. However, the operation fails, since the layer contains unsapported geometry types,
#' 
## ---- error=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x = st_interpolate_aw(stat, nafot, extensive = TRUE)

#' 
#' namely, `MULTISURFACE` geometries:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
as.data.frame(table(st_geometry_type(stat)))

#' 
#' We can "cast" the latter to `MULTIPOLYGON`, as follows:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
stat = st_cast(stat, "MULTIPOLYGON")

#' 
#' Now, all `r nrow(stat)` features in `stat` are of type `MULTIPOLYGON`:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
as.data.frame(table(st_geometry_type(stat)))

#' 
#' ## Fix topolopgy
#' 
#' It may also be necessary to fix the topology errors (depending on GEOS version):
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
stat = st_make_valid(stat)

#' 
#' ## Area-weighed interpolation
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x = st_interpolate_aw(stat, nafot, extensive = TRUE)
x$Group.1 = NULL

#' 
#' The result...
#' 
## ---- fig.width=9, fig.height=5, out.width="100%"--------------------------------------------------------------------------------------------------------------------------------------------------------
plot(x, max.plot = 16)

#' 
#' ## Calculating proportions
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dat = st_drop_geometry(x)
rownames(dat) = nafot$name_eng
dat[, 1:6]

#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dat = sweep(dat, 1, rowSums(dat), "/")
dat[, 1:6]

#' 
#' ## Hierarchical clustering
#' 
#' Now we can calculate the distance (i.e., the dissimilarity) between Nafot age structures, and apply hierarchical clustering:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
d = dist(dat)
hc = hclust(d, "average")

#' 
#' The dendrogram can then be split at a particular "level" to obtain a patricular number of distinct groups (such as `k=3`). The named vector `groups` specifies the classification of `nafot` features into three distinct clusters `1`, `2`, and `3`:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
k = 3
groups = cutree(hc, k = k)
groups

#' 
#' The dendrogram and the division to groups are shown in Figure \@ref(fig:hc-dendrogram):
#' 
## ----hc-dendrogram, fig.cap="Hierarchical dendrogram of age structure per Nafa", fig.width=7, fig.height=4, out.width="100%"-----------------------------------------------------------------------------
plot(hc)
rect.hclust(hc, k = k)

#' 
#' The division to clusters can also be displayed on a map (Figure \@ref(fig:hc-map)):
#' 
## ----hc-map, fig.cap="Age structure clusters on a map"---------------------------------------------------------------------------------------------------------------------------------------------------
nafot$group = factor(groups)
plot(nafot["group"])

#' 
#' Finally, to visualize the age structure distribution within each cluster, we can reshape the data:
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

#' 
#' and then use `ggplot2` (Figure \@ref(fig:ggplot-1)).
#' 
## ----ggplot-1, fig.cap="Distribution of age structures in each cluster of `nafot`", fig.width=8, fig.height=4, out.width="100%"--------------------------------------------------------------------------
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

#' 
#' # More information
#' 
#' [![half-size image](images/book1.png){width=45%}](https://geocompr.robinlovelace.net/)
#' [![half-size image](images/book2.png){width=45%}](https://www.r-spatial.org/book/)
#' [![half-size image](images/book3.png){width=45%}](https://geobgu.xyz/r/)
#' [![half-size image](images/book4.png){width=45%}](https://r-spatial.github.io/sf/articles/)
#' 
#' Other:
#' 
#' * `sf` tutorial from [useR!2017](https://edzer.github.io/UseR2017/) conference
#' * `sf` tutorial from [rstudio::conf 2018](https://edzer.github.io/rstudio_conf/#1) conference
#' * The r-spatial [blog](http://r-spatial.org/)
#' 
#' **Thank you for listening!**
