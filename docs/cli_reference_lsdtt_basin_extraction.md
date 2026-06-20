# lsdtt3 basin extraction guide 

This document explains the design choices of different basin objects in `lsdtt3`
and how to use these objects to extract basins.

---
## The `lsdtt_basin` object

The `lsdtt_basin` object holds a single basin, based off an outlet node defined
by a `flat_index` from `lsdtt_flow_base`, the flow info object. We assume that
the `flat_index` has already been checked to be valid.

The object holds an `lsdtt_flat_index_slice`, which is defined in 
`lsdtt_flow_types.hpp` and points to a slice in the flow info s vector (the 
Braun & Willett (2013) stack). The `lsdtt_basin` object therefore does not copy
any data. 

Note that the stack is the result of a depth-first search, which allows for this
slicing of upstream nodes. Barnes (2019, Geomorphology) argues that a breadth-
first search is better for parallel computing of an upstream-to-downstream D8
order, however, the resulting order does not allow for upstream node slicing.
Flow graphs resulting from multiple flow direction algorithms are also based
on the breadth-first search (Anand et al., 2020, Environmental Modelling & 
Software) and will therefore not work with the `lsdtt_basin`. 

The `lsdtt_basin` has an abstract interface, `lsdtt_basin_base`, which lives in
`src/core/`. This allows for a GPU version of `lsdtt_basin` to be developed in 
the future. The CPU version with the actual code lives in `src/cpu/` and is 
called `lsdtt_basin_cpu`. `lsdtt_basin_base` further contains shared types.

The `lsdtt_basin_cpu` will have various member functions:
- `check_basin_completeness(..)` takes boundary conditions and checks for each
pixel if it is on an open (nonperiodic, nonbaselevel) boundary and if it borders
a nodatavalue pixel. If this is true for any pixel, then the basin is incomplete
and the function returns false. 
- `get_basin_outline(..)` runs an edge tracing algorithm to get the outline of
the basin. It returns an lsdtt_basin_outline, which is yet to be designed.
- `get_basin_values(..)` takes an `lsdtt_raster_cpu` and uses the 
`lsdtt_flat_index_slice` to sample all basin values and put those in a vector,
which is returned. 
- `get_number_of_basin_pixels(..)` returns the number of pixels in the basin.

Functions that call a statistical method that is contained in `src/utils/`
should **not** be member functions of the `lsdtt_basin_cpu` object, but put them
instead in a different analysis file as free functions. This is to comply with
the design rule: "Write free functions unless you need private data members".

## The `lsdtt_junction_basin` object
The `lsdtt_junction_basin` is a child of `lsdtt_basin` that is linked to the 
`lsdtt_junction_network`. The `lsdtt_junction_network`, combined with the flow
info, defines the channel network. The `lsdtt_junction_basin` therefore has 
additional data members that denote the position of the basin in this network.

At any junction, two channels (or, exceptionally, up to seven, in a D8 flow 
routing world) combine that each define a basin. The pixels immediately upstream 
of that junction define the outlets of these basins. From these outlets, you can
go upstream and the first junction you encounter will be uniquely associated 
with that basin. So, aditional data members are an int32 recording the basin's
`junction_number` and the associated `pixel_stream_order`.

An `lsdtt_junction_basin` can still be selected by the outlet node, even if that
is not exactly a junction. The first junction encountered in upstream direction
is noted down and will correspond to the correct stream order.

## The `lsdtt_cosmo_basin` object 
The `lsdtt_cosmo_basin` is also a child of `lsdtt_basin`. It has additional data
members for holding cosmogenic information. It is used to process cosmogenic
nuclide samples, which are taken at one specific point that is not necessarily a
junction, so it does not inherit propeties from `lsdtt_junction_basin`.

## The `lsdtt_basin_manager` object (or structure)?
This is analysis code that lives in `/src/cpu`. It manages basin selection,
filtering and provides the opportunity to apply `/src/utils` statistical
analysis functions to all basins.

When making a basin manager, we need to decide whether we want to work with 
`lsdtt_basin` objects or with `lsdtt_junction_basin` objects. 

To select `lsdtt_basin`s, it needs to get a vector of `flat_index` values that
denote outlets. For each outlet, an `lsdtt_basin` is generated. In this case, 
the basin manager is only dependent on flowinfo.

To select `lsdtt_junction_basin`s, it needs to get an `lsdtt_junction_network`.
For each junction in the network, an `lsdtt_junction_basin` is generated. 
It should also be possible to get `lsdtt_junction_basin`s by supplying outlet
`flat_index` values. For each outlet, it should iterate up the flow info and
make `lsdtt_junction_basin` instances at each encountered junction. The
downstream-most basins (which initiated from a flat_index) should hold that
supplied `flat_index` as the outlet.

Then, there are a bunch of filtering options. These are member functions that
return `*this`, for optional chaining. Filtering methods are:
- `filter_incomplete(..)` calls the `check_basin_completeness` method of 
individual basins, and basins that return `false` are deleted.
- `filter_by_area(..)` takes a min and max number of pixels and removes those 
basins that do not have a number of pixels in this range.
- `filter_by_stream_order(..)` takes a min amd max stream order and removes 
those basins that are not in the range. Only works on an `lsdtt_junction_basin`
- `filter_nested(..)` removes all basins that are nested (upstream of another
basin)
- `filter_duplicates(..)` removes basins that have the same outlet `flat_index`
as another one.
- `filter_largest(..)` filters out all basins but the one with the largest
number of pixels.

For applying `/src/utils` functions, need to come up with something clever. 
All statistics functions expect a vector of data. `lsdtt_basin_cpu` has a
`get_basin_values(..)` option. So `lsdtt_basin_manager` can be used to pass
for each basin the value on to statistics functions. But I don't want to 
write an endless number of member functions to lsdtt_basin_manager, that 
just passes values into and out of the statistics functions. need to come up
with a clever way for a cli to use the basin manager and appy the statistics
of choice to it, and writes it to a file as desired by the cli using the
lsdtt_data_loader.

## Selecting basins that cross a line (e.g. a fault)

`select_basins_from_line_crossings` is a basin-selection mode for studying basin
geometry relative to a baseline structure such as a fault or a valley axis. It
finds the channel pixels that lie on a supplied polyline and whose contributing
area exceeds a threshold, treats each as a basin outlet, and extracts the basin
upstream of it. Because the analysis is stream-order / geometry based, it
automatically keeps only the **largest complete basin per independent crossing**
(it forces the incomplete- and nested-basin filters on).

The supplied line need not be in the DEM CRS — it is reprojected into the DEM
CRS internally (its layer CRS is read and a coordinate transform applied), so
you can pass, for example, a WGS84 fault trace straight from the GEM database or
from `lsdtt-raster-preprocessing`'s `clip_vector_data_to_raster` output.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `select_basins_from_line_crossings` | bool | `false` | Enable the mode. Mutually exclusive with the other `select_*` modes. |
| `crossing_line_fname` | string | `NULL` | Polyline vector the channels must cross (any OGR format; any CRS). Required. |
| `crossing_minimum_contributing_pixels` | int | `500` | Minimum contributing pixels at a crossing to seed a basin outlet. |

Combine with `filter_basins_by_area` + `basins_max_area_pixels` to drop a trunk
valley that the line runs *along* (rather than across), keeping only the
fault-transverse catchments.

### Per-basin analysis output

In this mode the tool writes `<write_fname>_line_crossing_basins.csv` — one row
per retained basin, with every per-basin quantity computed in C++ (outlet, centre
of mass, area, the outlet→COM bearing, and the bearing relative to the local
fault orientation):

| Column | Description |
|---|---|
| `basin_id` | Index of the basin in the retained set |
| `area_m2` | Basin area (upslope of the crossing), m² |
| `outlet_x`, `outlet_y` | Crossing-pixel (outlet) cell-centre, raster CRS |
| `com_x`, `com_y` | Basin centre of mass, raster CRS |
| `outlet_to_com_bearing_deg` | Bearing from outlet to COM (deg from North, clockwise) |
| `local_fault_azimuth_deg` | Azimuth of the nearest line segment to the basin |
| `side` | Which side of the line the basin lies on (`east`/`west`, by COM easting) |
| `deflection_from_fault_orthogonal_deg` | Outlet→COM bearing minus the line-orthogonal pointing into the basin, in (−180, 180]. ~0 means the basin extends straight across the line; a systematic sign by side indicates shear-related bending. |

The polygons (`write_basin_outlines_as_polygon`) and this CSV are consumed by the
`lsdviztools3` `fault_basin_deformation.py` example for plotting (rose diagrams of
basin bearings, deflection-by-side distributions, etc.).

Example parameter file:

```
read_path:   /data/parkfield/
write_path:  /data/parkfield/
read_fname:  parkfield_COP30_fill.tif
write_fname: parkfield
raster_is_filled: true
source_pixels_threshold: 200

select_basins_from_line_crossings: true
crossing_line_fname: parkfield_COP30_clipped_vector.gpkg
crossing_minimum_contributing_pixels: 500

# focus on fault-transverse catchments (drop the along-fault trunk valley)
filter_basins_by_area: true
basins_max_area_pixels: 22222

write_basin_outlines_as_polygon: true
```
