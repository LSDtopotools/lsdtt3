"""Compact, retrieval-oriented API context for language models.

This module intentionally favours exact signatures and behavioural constraints
over prose.  It is not human-facing reference documentation: it is what an AI
assistant reads when it has only the installed package to go on (for example
a student's Colab session), so every name and argument in it is checked
against the public API manifest by ``tools/check_llm_cards.py``.

Authors
-------
B. Gailleton
S. M. Mudd
"""

from __future__ import annotations

from collections.abc import Mapping
import inspect
import re

from ._llm_status import EXPERIMENTAL_OVERRIDES, READY
from ._topics import GUIDE_TOPICS, guide_topic


class _PlainText(str):
    """A ``str`` that displays as its text rather than as an escaped literal.

    ``helper_for_llm()`` is usually the last expression in a notebook cell,
    where a plain ``str`` is shown as one long quoted line full of ``\\n``.
    This subclass is still a ``str`` in every other respect.

    Authors
    -------
    S. M. Mudd
    """

    __slots__ = ()

    def __repr__(self) -> str:
        return str(self)

    def _repr_pretty_(self, printer, cycle) -> None:
        # IPython (Jupyter, Colab) asks for this before falling back to repr().
        printer.text(str(self))


# The core workflow. The module docstring of lsdtt3 carries the same code;
# tests/test_llm_helper.py checks that the two copies agree.
WORKFLOW_CODE = '''\
import lsdtt3

# 1. Download a 30 m Copernicus DEM for a WGS84 box (west, south, east, north).
#    It comes back projected to the UTM zone of the box centre, in metres.
dem = lsdtt3.fetch_dem(15.63, 37.88, 16.21, 38.31, source="cop30_aws", grid_spacing=30.0)
dem = lsdtt3.mask_to_nodata(dem, 0.0)  # Copernicus stores the sea as 0 m
filled = lsdtt3.fill_pits(dem)

# 2. D8 flow routing and the channel network.
threshold = 500  # contributing pixels needed to start a channel
F = lsdtt3.FlowInfo(filled)
N = F.channel_network(F.channel_scope(threshold_contributing_pixels=threshold))

# 3. Snap (latitude, longitude) outlets onto the channel network: each point
#    follows the flow downstream to the first pixel that drains at least
#    `threshold` pixels, as basin_outlet_snapping_method: downstream does in
#    the CLIs. The result is an N x 2 array of (row, col).
outlet_pixels = F.snap_points(
    latitude=[37.94234, 38.06143],
    longitude=[15.79008, 16.06060],
    method="downstream",
    threshold_contributing_pixels=threshold,
)

# 4. Basins: one per outlet, basin_id 0, 1, ... in outlet order.
B = F.basins(outlet_pixels)
print(B.table[["basin_id", "outlet_row", "outlet_col", "area_metres2", "complete"]])
polygons = B.polygons  # GeoDataFrame in the DEM's CRS

# 5. Chi and normalised channel steepness (ksn) by segment fitting.
fit = N.chi(outlets=outlet_pixels, m_over_n=0.45).fit(filled, critical_divergence_metres=10.0)
print(fit.segments[["basin_id", "segment_id", "ksn"]])
print(fit.summary[["basin_id", "median_ksn"]])

# 6. Nothing is written to disk unless you ask.
filled.write("filled.tif")
lsdtt3.hillshade(filled).write("hillshade.tif")
polygons.to_file("basins.gpkg")
'''


_INDEX = """LSDTT3/API3 compact lsdtt3 guide for AI assistants (import lsdtt3; names below are lsdtt3.<name>)
How to work: 1. read helper_for_llm('gotchas'); 2. read the topic for the task; 3. before calling anything, get its exact arguments with helper_for_llm('symbol:<name>') or 'symbol:<Class>.<method>'; 4. if no topic fits, helper_for_llm('search:<words>') (a CLI parameter name works too). Topics list the ready exports; symbol and search also find the others, marked [experimental].
Notation: R=Raster M=RasterMeta F=FlowInfo G=FlowGraph H=RasterHierarchy S=ChannelScope N=ChannelNetwork B=BasinCollection C=ChiCollection L=LandscapeModel Q=CRNColumn V=geopandas.GeoDataFrame D=pandas.DataFrame; name?=optional argument; ...=more optional keywords; ->=returns; a|b=either.
Rules: projected grids in metres; slope=gradient m/m; A=R.to_numpy() is A[row,col] with row 0 at the NORTH edge; R.mask_numpy() True=valid; results are new in-memory objects, save with R.write(path); V must use R.metadata.crs_wkt; outlets are (row,col) arrays: get them from lat/lon or x/y with F.snap_points(...).
Topics:
  workflow: runnable fetch DEM -> fill -> flow -> channels -> basins -> chi/ksn -> write
  gotchas: the common mistakes (sea stored as 0 m, CRS, row/col, basin ids, column names)
  terrain: rasters and grid metadata, conditioning (fill, breach), surface metrics, hillshade
  fetch: DEMs and bounds from the web, UTM zones
  flow: FlowInfo (D8 routing), drainage area, snapping points (FlowGraph: topic drainage)
  channel: channel scope and network, channel heads, slope-area and low concavity
  basin: basins, divides, basin metrics
  chi: chi, ksn, segment fits, knickpoints, inversion
  vector: GeoDataFrames: rasterize, clip, sample, swath profiles, valley centrelines and sections
  water: river lines, burning lines into a DEM, water centrelines, line repair
  hydraulic: GraphFlood
  landscape: landscape evolution, hillslope transport
  sediment: sediment routing, landslides
  cosmo: cosmogenic nuclides
  drainage: FlowGraph (D-infinity, MFD area), depressions, internal drainage, wetness and connectivity indices
  synthetic: synthetic landscapes and noise
  stats: raster arithmetic, thresholds, comparisons
Also: 'all' = this index, the guidance of every topic and the name of every ready export, in one text for pasting into a chat; 'catalogue' = one line per ready export.
"""

_WORKFLOW = (
    "@workflow fetch DEM -> remove sea -> fill -> flow -> channels -> snap outlets"
    " -> basins -> chi/ksn -> write. Runnable as is (fetch_dem needs network).\n"
    "```python\n" + WORKFLOW_CODE + "```\n"
    "B.table: basin_id outlet_flat outlet_row outlet_col n_pixels area_metres2 complete"
    " outlet_on_dem_edge parent_id nesting_depth (complete=False: basin cut by the DEM edge"
    " or nodata, so enlarge the box).\n"
    "fit.nodes: basin_id outlet_flat flat receiver_flat chi observed_elevation fitted_elevation"
    " residual_elevation branch_id is_main_stem is_source is_junction source_key (one row per channel pixel).\n"
    "fit.segments: basin_id segment_id parent_segment_id branch_id downstream_chi upstream_chi"
    " downstream_elevation upstream_elevation ksn n_nodes source_key (one row per straight chi segment).\n"
    "fit.summary: basin_id n_segments mean_ksn median_ksn min_ksn max_ksn chi_weighted_mean_ksn"
    " main_stem_mean_ksn main_stem_median_ksn (one row per basin).\n"
)

_GOTCHAS = """@gotchas
- Copernicus DEMs (fetch_dem source='cop30_aws') store the sea as 0 m, not nodata. Run mask_to_nodata(R,0.0) (sets values <= 0 to nodata) before fill_pits and FlowInfo, or the sea becomes one huge flat basin.
- Results are new in-memory objects; inputs are never changed and nothing is written unless you call R.write(path) (or V.to_file(path) for GeoDataFrames).
- Install with pip install lsdtt3 (Colab: !pip install lsdtt3); it needs Python 3.12 or later. There is one wheel per platform (cp312-abi3 in the file name); if you download one from a GitHub release, keep its file name, pip reads the Python and platform from it. The wheel also installs the lsdtt-* command-line programs (e.g. lsdtt-basin-extraction), which are driven by .param files.
- The 0.7.1 wheel sets PROJ_DATA, PROJ_LIB and GDAL_DATA when imported (later versions do not). If pyproj or GMT (pygmt, lsdviztools3) then fail, remove them straight after import lsdtt3: for name in ('PROJ_DATA','PROJ_LIB','GDAL_DATA'): os.environ.pop(name,None)
- Grids must be projected with metre units; fetch_dem already reprojects to UTM. Use reproject(R,epsg) for a geographic raster. Slopes are gradients (m/m), never degrees.
- A=R.to_numpy() is A[row,col] and row 0 is the NORTH edge. Pixel centre: x=M.x_min+(col+0.5)*M.cell_size, y=M.y_max-(row+0.5)*M.cell_size with M=R.metadata; M.row_col_to_xy(row,col) and M.xy_to_row_col(x,y) do this for arrays (-1 off the grid), and M.lonlat_to_row_col(longitude,latitude) takes WGS84 degrees. Invalid pixels hold the nodata value: use R.mask_numpy().
- GeoDataFrame inputs (polygons, lines, points, swath baselines) must already be in the raster's CRS: V=V.to_crs(R.metadata.crs_wkt).
- F.basins(outlets) takes an N x 2 integer array of (row, col) pixels, not x/y or lat/lon, and silently drops duplicate outlets, so B.size can be smaller than N. Convert and snap points first with F.snap_points(latitude=...,longitude=...) or F.snap_points(x,y) (see workflow); do not hand-roll the row/col arithmetic. snap_points methods and thresholds are the CLIs' basin_outlet_snapping_method and basin_outlet_contributing_pixels_threshold: pass the pixel count as threshold_contributing_pixels, or an area in m^2 as threshold_drainage_area_m2 (not both); with neither the threshold is the CLI default, 1000 contributing pixels, for every method. snap_points raises ValueError naming every point that is off the grid or cannot be snapped (with the C++ snap status); widen max_search_radius_pixels, move the point or drop it.
- Basin IDs are zero-based in tables (B.table, B.polygons, fit.nodes/segments/summary), but the raster B.labels stores basin_id+1, with 0 outside every basin. Use B.mask(basin_id) for one basin.
- fit.nodes has a 'flat' pixel index but no row/col: row=flat//n_cols, col=flat%n_cols, or F.row_col(flat_index). The same holds for every *_flat column.
- F.drainage_area() is in m^2; F.contributing_pixels() is a pixel count, as are thresholds such as threshold_contributing_pixels.
- Column names differ from the lsdtt-* CLI CSVs that lsdviztools3 expects: rename basin_id->basin_index, and ksn->k_sn in fit.segments, before using lsdviztools3 chi plots: fit.segments.rename(columns={'basin_id':'basin_index','ksn':'k_sn'}). fit.tributaries and fit.knickpoints need only basin_id->basin_index.
"""

# The hand-written part of each topic. A catalogue topic (GUIDE_TOPICS) is
# followed by its generated catalogue when rendered (_topic_card).
_GUIDANCE = {
    "workflow": _WORKFLOW,
    "gotchas": _GOTCHAS,
    "terrain": """@terrain rasters, conditioning, surface metrics
Raster.read(path)->R; Raster.from_numpy(data,mask,metadata:M)->R (data float32 [row,col]; mask bool, True=valid)
R.write(path,format='cog'); R.to_numpy()->float32 array copy; R.mask_numpy()->bool array; R.metadata->M; R.shape->(n_rows,n_cols); R.valid_pixel_count
RasterMeta.from_shape(n_rows,n_cols,x_min=0,y_min=0,cell_size=1,nodata_value=-9999,crs_wkt='',epsg_code=0)->M; M.x_min M.y_max M.cell_size M.crs_wkt M.epsg_code M.n_rows M.n_cols
M.xy_to_row_col(x,y)->(row,col) (-1 off the grid); M.row_col_to_xy(row,col)->(x,y) pixel centres; M.lonlat_to_row_col(longitude,latitude,crs='EPSG:4326')->(row,col); all take scalars or arrays
mask_to_nodata(raster:R,threshold,below=True)->R; fill_pits(raster:R,min_slope=0.0001,...)->R; breach_depressions(raster:R)->R; breach_least_cost(raster:R,...)->R
hillshade(raster:R,azimuth=315,altitude=45,z_factor=1)->R (0-255); multidirectional_hillshade(raster:R,n_directions=8,...)->R; slope(raster:R)->R (3x3 gradient, m/m)
polyfit_metrics(raster:R,metrics,lambda_metres=0,...)->dict of R keyed by metric; metrics is a list such as ['slope','aspect_degrees','plan_curvature','profile_curvature','tangential_curvature','mean_curvature','maximum_curvature','minimum_curvature','gaussian_curvature','laplacian','dzdx','dzdy','fitted_surface']; 'slope' is a gradient (m/m); lambda_metres is the fit window (0 = three cells)
classify_surface_shape(raster:R,lambda_metres=0,...)->dict; classify_surface_shape(maximum_curvature:R,minimum_curvature:R,kt=0)->dict
reproject(raster:R,epsg,cell_size=0,resampling='bilinear')->R; reproject_to_grid(raster:R,target_metadata:M)->R; mosaic(rasters)->R; clip_to_polygons(raster:R,polygons:V)->R
resample_raster(source:R,target_metadata:M,method='linear')->R; double_raster|halve_raster(source:R,method)->R; pad_to_nearest_power_of_two(source:R)->R
H=RasterHierarchy(finest:R,storage='retain_allocated',keep_finest=True,levels=0); H.restrict_to(level,method='conservative_mean')->R; H.prolong_to(level,method='linear')->R
spectral_analysis(raster:R,...)->SpectralResult; frangi_filter(raster:R,...)->FrangiResult; swath_profile(reference:R,baseline:V,values=R?,half_width_metres=0,bin_width_metres=1000)->SwathProfile
""",
    "fetch": """@fetch download DEMs (needs network); results are projected in-memory rasters
fetch_dem(west,south,east,north,source='cop30_aws',target_epsg=0,grid_spacing=30.0,...)->R
  west,south,east,north are WGS84 decimal degrees (four floats; not a bounds object or a GeoDataFrame). source: 'cop30_aws' (keyless Copernicus 30 m) | 'opentopography' (needs OT_API_KEY) | 'planetary_computer'. target_epsg=0 means the UTM zone of the box centre. grid_spacing is in metres.
b=bounds_from_vector(path,layer='')->VectorBounds (path is a vector FILE with polygons or exactly two corner points); b.ok; b.bbox is the tuple (west,south,east,north), so use fetch_dem(*b.bbox)
fetch_dem_for_outlets(outlets:V,source='cop30_aws',grid_spacing=30.0,...)->OutletFetchResult (.rasters list of R, .groups, .snapped_outlets, .basin_labels)
utm_epsg_for_point(longitude,latitude)->int; cop30_tiles_for_bbox(west,south,east,north)->list[str]; opentopography_datasets()->list; http_available()->bool
Oblique frames (a grid aligned with a fault/river/range, feature along +x, so the strip is not padded to its bounding box):
  frame=RotationFrame.from_feature(path,width=0.0,end_buffer=-1.0,grid_spacing=0.0,layer='',source=None) (path: line, two points or polygon FILE; width in m for lines)
  fetch_dem(frame=frame,grid_spacing=frame.cell_size)->R in the frame; frame.write(path); RotationFrame.read(path)
  existing DEM: frame=RotationFrame.from_feature(path,width=1200.0,source=dem); rot=rotate_to_frame(raster=dem,frame=frame)->R; or RotationFrame.from_valid_data(source=dem)
  back to the source grid: rotate_back(raster=basins,frame=frame,resampling='near')->R (rotate back products, never the DEM)
  bearings (aspect, link bearings, ...) are from GRID north (90 = along the feature) unless bearings_relative_to_true_north=True
! Copernicus stores the sea as 0 m: dem=mask_to_nodata(R,0.0). Save explicitly: R.write(path).
""",
    "flow": """@flow D8 routing (FlowInfo); multiple-direction routing (FlowGraph) is in topic drainage
F=FlowInfo(dem:R,...): pass a pit-filled DEM (fill_pits). F is immutable D8 topology; boundary_conditions has one character per edge, north east south west.
Units: F.drainage_area()->R (m^2); F.contributing_pixels()->R (pixel count); F.flow_distance_to_outlet()->R (metres). Pixels are (row,col) or flat=row*n_cols+col (F.row_col, F.flat_index).
F.receiver_of(row,col)->(row,col) or None at an outlet
F.snap_points(x?,y?,longitude?,latitude?,points?,crs?,method='downstream'|'expanding_square'|'largest_area',threshold_contributing_pixels?,threshold_drainage_area_m2?,...)->N x 2 int64 (row,col) array for F.basins; coordinates are keyword-only except x,y; same rules as the CLIs' basin_outlet_snapping_method: 'downstream' = follow D8 receivers to the first pixel with contributing pixels >= threshold, 'expanding_square' = nearest such pixel in rings min..max radius, 'largest_area' = largest drainage area within max_search_radius_pixels; threshold default 1000 contributing pixels (CLI default); threshold_drainage_area_m2 (m^2) becomes the smallest N pixels with N*pixel_area >= it; raises ValueError for points off the grid or not snapped
Chain: S=F.channel_scope(threshold_contributing_pixels=500,...) or F.scope_from_heads(heads); N=F.channel_network(S); F.basins(outlets?,min_upslope_pixels=0)->B (outlets: N x 2 int array of (row,col); None = every D8 outlet)
F.height_above_channel(dem:R,channels:R)->R (HAND; channels is a 0/1 raster such as N.channel_mask)
""",
    "channel": """@channel channel heads and networks
S=F.channel_scope(threshold_contributing_pixels=...), or heads=extract_channel_heads(dem:R,flow:F,method='threshold'|'wiener',threshold_contributing_pixels=1000,...)->ChannelHeads then S=F.scope_from_heads(heads)
N=F.channel_network(S); N.channel_mask->R (0/1); N.n_channel_cells; N.n_junctions
N.lines_geodataframe()->V (one line per link, with stream_order); N.points_geodataframe()->V; N.junctions_geodataframe()->V; N.links()->D; N.junctions()->D; N.junction_angles()->D; N.link_profiles(values:R)->D
C=N.chi(outlets?,m_over_n=0.45,reference_area=1.0,denominator?); N.chi_sweep(dem:R,m_over_n_values,...)->D
F.channel_point_attributes(filled_dem:R,points:V,*,threshold_contributing_pixels,slope_window_half_length_metres=25.0,channel_points_snap_max_search_radius_pixels=10)->V (CLI run_channel_point_attributes: each point moved to the NEAREST channel cell by straight-line distance, not along D8; keeps every column, adds original_x,original_y,snap_distance_m,snap_status,snapped_flat,contributing_pixels,drainage_area_m2,flow_distance_m,elevation_m,smoothed_elevation_m,slope_m_per_m (along-channel least-squares slope over +/- the half-length),slope_valid,slope_window_length_m; geometry = snapped cell centre in the grid CRS; unsnapped points keep position with NaN; a clashing added name gets _lsdtt)
find_channel_candidates(raster:R,...)->ChannelCandidatesResult
Slope-area (CLI lsdtt-slope-area; same routines, same defaults): carve then fill, and leave the changed pixels out: routing=fill_pits(breach_least_cost(raster:R,fill_remainder=False)); mask=conditioned_pixel_mask(input_dem:R,routing_dem:R) (0 untouched, 1 filled, 2 carved; pass it as exclude=mask); F=FlowInfo(dem:R) and area=FlowGraph(dem:R,method='dinf').drainage_area(), both on routing (D-infinity, the CLI default) or F.drainage_area() (D8)
gradient: polyfit_metrics(raster:R,metrics,lambda_metres=150)['slope'] on the INPUT dem (hillslope-channel transition), or F.d8_flow_path_gradient(routing_dem:R,drainage_area:R,window_metres=150,exclude=R) (along channels; needs window >= 4 cells); slope_area_bins(drainage_area:R,slope:R,n_log_bins=20,exclude=R)->D (the CLI's lambda CSV columns)
Low concavity: F.low_concavity(drainage_area:R,gradient:R,threshold_contributing_pixels=1000,max_abs_concavity=0.15,exclude=R)->D (one row per channel pixel; low_concavity 1, 0, or -1 for no valid fit; gradient must be the D8 flow-path gradient); F.low_concavity_segments(...,min_run_pixels=1)->V. A flag means S did not vary with A in the window, not a clean power law: check theta_se and r2. With D8 area, flags gather at confluences.
""",
    "basin": """@basin drainage basins
B=F.basins(outlets,min_upslope_pixels=0) (outlets: N x 2 int array of (row,col), e.g. from F.snap_points(latitude=...,longitude=...); duplicate outlets are dropped)
B.table->D (basin_id from 0, outlet_row, outlet_col, n_pixels, area_metres2, complete, parent_id, ...); B.polygons->V (basin_id, geometry; DEM CRS); B.size
B.labels->R (basin_id+1; 0 outside every basin); B.mask(basin_id)->R; B.pixels(basin_id); B.membership(row,col)->basin_ids; B.statistics(dem:R,...)->D
B.outline_ridgelines(dem:R,split_method='longest_channel')->V (left and right ridge lines of each basin outline)
junction_basins(network:N,selection='all',outlets?)->JunctionBasinCollection (.table .labels .filter_area(min_pixels,max_pixels) .filter_incomplete() .filter_nested())
line_crossing_basins(network:N,line:V,minimum_contributing_pixels=500)->B; side_tributary_basins(flow:F,points_path,...)->SideTributaryBasins
basin_surface_metrics(basins:B,dem:R,scales,metrics,...)->dict; basin_channel_lengths(basins:B,network:N,max_stream_order=6)->D; basin_ridgeline_extrema(basins:B,network:N,dem:R,...)->V
divide_migration(flow:F,dem:R,outlets,...)->DivideMigrationResult
""",
    "chi": """@chi chi and normalised channel steepness (ksn)
C=N.chi(outlets?,m_over_n=0.45,reference_area=1.0,denominator?) (outlets: the same (row,col) array given to F.basins; None = every network outlet); C.n_basins; C.nodes()->D; C.raster()->R (chi)
fit=C.fit(dem:R,critical_divergence_metres=5.0,mainstem='max_flow_length',knickpoint_threshold_algorithm_choice='p95',run_step_change_detection=False,step_change_min_gradient=0.40,step_change_max_window_length_metres=90.0,step_change_min_drop_metres=5.0)->ChiFit (pass the filled DEM; 10 m suits 30 m data); ksn is the slope of each straight segment in chi-elevation space
fit.nodes->D (one row per channel pixel: basin_id, flat, chi, observed_elevation, fitted_elevation, branch_id, is_main_stem; no row/col, use F.row_col(flat_index))
fit.segments->D (one row per segment: basin_id, segment_id, ksn, downstream_chi, upstream_chi, n_nodes, source_key); fit.summary->D (one row per basin: median_ksn, mean_ksn, chi_weighted_mean_ksn, main_stem_median_ksn)
fit.tributaries->D (one row per tributary = fit branch: source_key, junction_drainage_area_m2, area_fraction, area_rank_in_basin, ksn_max, n_descendants); fit.knickpoints->D (one row per segment boundary: delta_ksn = up - down, at_junction, x, y, chi, elevation, split_delta_ksn_normalised = ranking measure, is_major); fit.knickzones->D (foot/lip pairs; foot_type split|basin_outlet|tributary_junction: foot_chi, lip_chi, relief_m, ksn_inside); fit.steps->D (run_step_change_detection: waterfall/steep reaches: bottom_flat, top_flat, drop_m, at_tributary_mouth, at_basin_outlet); fit.step_aware_fit->ChiFit|None (stage 2, cut at steps: k_sn only, nodes/segments gain is_step; knickpoint tables raise, rank on stage 1; stage 1 knickpoints gain nearest_step_id, nearest_step_chi_distance, step_in_adjacent_segments); fit.knickpoint_thresholds->D (p90/p95/median_plus_3mad/otsu per |measure|; q3_plus_3iqr on ksn)
fit.disorder(); fit.disorder_bootstrap(replicates=1000,...); inv=fit.invert_ksn(m_over_n_values,...)->ChiInversion; inv.nodes; inv.summary; inv.u_star_raster(basin_id=0)->R
N.chi_sweep(dem:R,m_over_n_values,outlets?,...)->D (disorder per basin for each m/n); F.chi_from_outlets(outlets?,m_over_n=0.45,...)->R (chi raster without a network)
! chi needs N built from F on the same grid as the elevation R. With reference_area=1, ksn has units of m^(2*m_over_n).
""",
    "water": """@water mapped water lines: fetching river lines, burning them into a DEM, water centrelines, line repair; V must be in the raster CRS
burn_streams(dem:R,lines:V)->StreamBurnResult (.raster, .profiles); fetch_river_lines(river_line_source='gb_os_open_rivers'|'local_file',west,south,east,north or bounds_fname=dem path,...)->RiverLinesResult (.lines downstream-oriented V, .nodes, .topology, .attribution; CLI lsdtt-fetch-river-lines names); burn_to_corridor(dem:R (not filled),lines:V,*,vector_name_field='',vector_corridor_half_width_metres=0,burn_profile_max_pit_pixels=16,burn_profile_max_cut_metres=2000,burn_profile_fit_algorithm_choice='linear'|'chi',...,min_elevation=-9999)->BurnCorridorResult (.dem final conditioned R, .burned, .burn_depth, .profiles D, .channels/.exit_routes/.cross_path_receivers V, .report with downhill checks); vector is the truth: profile through the DEM minima under each line, rest of DEM carved+filled to drain to it; CLI lsdtt-stream-burning stream_burn_mode burn_to_corridor; water_centrelines(areas:V,tributaries:V,...)->WaterCentrelineResult; repair_water_lines(water_lines:V,roads:V,...)->WaterLineRepairResult
""",
    "vector": """@vector GeoDataFrame inputs and outputs; V must be in the raster CRS: V=V.to_crs(R.metadata.crs_wkt) (water lines and burning: topic water)
rasterize_polygons(raster:R,polygons:V,coverage='majority')->R (0/1); clip_to_polygons(raster:R,polygons:V,coverage='majority',crop=True)->R
sample_points(raster:R,points:V,column='raster_value')->V; linear_stats(raster:R,lines:V,spacing_metres=10)->V; mask_from_points(raster:R,points)->R
extract_valley_centrelines(dem:R,scope:S|channels,...,valley_mask=R)->ValleyCentrelineResult (valley_mask is keyword-only and required); extract_valley_cross_sections(dem:R,centrelines,...)->ValleyCrossSections
swath_profile(reference:R,baseline:V,...)->SwathProfile (baseline is exactly one LineString)
bounds_from_vector(path,layer='')->VectorBounds (reads a vector FILE, not a GeoDataFrame)
""",
    "hydraulic": """@hydraulic
fill_hydraulic_surface(surface:R,node_state:R,boundary_conditions='nnnn',...)->HydraulicSurfaceFill
run_graphflood(bed:R,node_state:R,rainfall_rate:R,friction_coefficient:R,...)->HydraulicResult
run_graphflood_optimised(bed:R,node_state:R,rainfall_rate:R,friction_coefficient:R,...)->GraphFloodOptimisedResult
run_graphflood_optimised_v2(bed:R,node_state:R,rainfall_rate:R,friction_coefficient:R,prescribed_qin:R|None=None,...)->GraphFloodOptimisedResult
run_graphflood_multigrid(bed:R,node_state:R,rainfall_rate:R,friction_coefficient:R,backend='optimised_v2',...)->GraphFloodMultigridResult
run_graphflood_induced(bed:R,node_state:R,rainfall_rate:R,friction_coefficient:R,qin_threshold,...)->GraphFloodInducedResult
node_state: 0=outside/nodata,1=active,2=flux-only ghost,3=outlet; ghosts route prescribed_qin but evolve depth/qout from initial depth, only 3 drains; each active component needs an outlet unless check_valid_outlet=False.
! hydraulic routing is on z+h; stationary workflow fills h, not bed elevation.
""",
    "landscape": """@landscape landscape evolution and hillslopes (sediment routing and landslides: topic sediment)
L=LandscapeModel(dem:R,incision_threshold_specific_area=0,fluvial_method='d8'|'dinf'|'none',hillslope_model='none'|'linear'|'roering'|'roering_implicit'|'cidre'|'adi'|'cellular',coupling_regime='automatic'|'coupled'|'channels_as_boundaries',...)
L.step(compute_erosion_rate=False)->LandscapeStepResult; L.run(n_steps,compute_erosion_rate=False)->LandscapeStepResult; L.set_forcing(uplift_rate,baselevel_fall_rate=0); L.elevation->R; L.time_years
erodibility>0 (K=0 raises; use fluvial_method='none'); erodibility_field=R|ndarray(n_rows,n_cols) absolute K>0 on the DEM grid; L.set_erodibility_field(R|ndarray|None) next step; L.erodibility_field->R|None
hillslope_transport_implicit(dem:R,channels:R,duration_years,diffusivity=0.001,uplift_rate=0,critical_slope=0.8,nonlinear=True,...)->HillslopeTransportResult (channels: 0/1 mask of fixed pixels)
hillslope_transport_explicit(dem:R,channels:R,diffusivity,duration_years,nonlinear=False,...)->R; hillslope_transport_adi(dem:R,diffusivity,duration_years)->R; hillslope_transport_cellular(dem:R,channels:R,duration_years,...)->HillslopeCellularResult
hillslope_steady_state(dem:R,channels:R,model='roering',...)->HillslopeSteadyStateResult
""",
    "sediment": """@sediment sediment routing and landslides
route_sediment(network:N,dem:R,formula='wilcock_crowe',...)->SedimentRoutingResult; trigger_landslides(dem:R,flow:F,target_volume_m3,...)->LandslideResult; run_landslide_series(dem:R,flow:F,...)->LandslideSeriesResult
""",
    "cosmo": """@cosmo cosmogenic nuclides (Be10, Al26, C14)
Q=CRNColumn(rock_density_kg_m3=2650); Q.set_erosion_history(time_yr_bp,erosion_rate_m_yr,removal_metres); Q.transient(sampling_depths_metres,...)->list[CRNColumnSample]; Q.steady_state(...)->CRNColumnSample
crn_production_scaling(dem:R,crn_parameters_path='',scaling='stone',nuclide='Be10')->R; crn_uniform_raster(R,value)->R
crn_steady_state(erosion_rate:R,production_scaling:R,nuclide='Be10',...)->R (atoms/g); crn_step_change(time_since_change:R,production_scaling:R,...)->R
crn_transient_history(column,production_scaling:R,topographic_shielding:R?)->{'be10':R,'al26':R,'c14':R}
crn_invert_basin_erosion(flow:F,samples:V,production_scaling:R,concentration_column='concentration',nuclide='Be10',...)->D; crn_outlet_concentration(concentration:R,erosion_rate:R,flow:F,row,col)->float; crn_pair_diagnostics(samples,...)
""",
    "drainage": """@drainage multiple-direction routing (FlowGraph), depressions, internal drainage, and wetness and connectivity indices
G=FlowGraph(dem:R,method='mfd'|'dinf',...): G.drainage_area()->R spreads flow over hillslopes
topographic_wetness_index(drainage_area:R,slope:R)->R and soil_topographic_index(...), sediment_connectivity_index(flow:F,roughness:R,slope:R,drainage_area:R,area_threshold_m2) (slope is a gradient, m/m)
find_internal_drainage(dem:R,method='hierarchy',...)->InternalDrainageResult (.candidates D, .polygons V); depression_hierarchy(dem:R,boundary_conditions='oooo',...)->DepressionHierarchyResult (.labels R, .flow_direction R, .nodes D)
fill_pits_boundary_aware(raster:R,boundary_conditions='oooo',...)->R; fill_pits_with_outlets(raster:R,outlet_mask:R,...)->R; fill_pits_with_sinks(raster:R,sinks:R|array,...)->R; make_outlet_mask(raster:R,north=False,...)->R
""",
    "synthetic": """@synthetic synthetic landscapes on a RasterMeta grid
synthetic_noise(metadata:M,kind='red',seed=1,...)->R; synthetic_mountain_range(metadata:M,...)->dict; synthetic_hillslope_strip(metadata:M,...)->dict
synthetic_network(metadata:M,outlet_mode='single_south',...)->SyntheticNetworkResult (.elevation R); add_uniform_noise(raster:R,seed,magnitude)->R
""",
    "stats": """@stats statistics and raster arithmetic
linear_stats(raster:R,lines:V,spacing_metres=10)->V; compare_bearings(first,second,axial=False,...)->dict; relief_ratio(raster:R,...)->R
add|subtract|multiply|divide(a:R,b:R)->R (two rasters on the same grid; no scalars); scale(raster:R,factor)->R (multiply by a number); absolute(raster:R)->R
threshold_above|threshold_below(raster:R,threshold)->R (0/1 mask); mask_by_binary(values:R,keep_mask:R)->R; invert_binary_mask(raster:R)->R
B.statistics(dem:R,...)->D; sample_points(raster:R,points:V)->V
""",
}

# Order of topics in helper_for_llm("all"): the workflow and pitfalls first.
_ALL_ORDER = ("workflow", "gotchas") + GUIDE_TOPICS

#: Classes whose method names ``all`` lists, since most work is done through
#: their methods (Simon, 2026-11-19).
_ALL_KEY_CLASSES = ("FlowInfo", "ChannelNetwork", "BasinCollection", "ChiCollection", "Raster")

_ALIASES = {
    "raster": "terrain", "surface": "terrain", "channels": "channel",
    "networks": "channel", "basins": "basin", "vectors": "vector",
    "graphflood": "hydraulic", "hydraulics": "hydraulic", "lem": "landscape",
    "landslides": "sediment", "landslide": "sediment", "burning": "water", "river_lines": "water", "water_lines": "water", "crn": "cosmo", "cosmogenic": "cosmo",
    "depression": "drainage", "internal": "drainage", "noise": "synthetic",
    "example": "workflow", "quickstart": "workflow", "pitfalls": "gotchas",
    "ksn": "chi", "fetching": "fetch", "statistics": "stats",
}


def _guide_status(name: str, namespace: Mapping[str, object]) -> str:
    """Return ``"ready"`` or ``"experimental"`` for an export or ``"Class.method"``.

    An export is ready if it is in ``READY`` and not in
    ``EXPERIMENTAL_OVERRIDES``; a method takes its class's status unless it is
    demoted on its own. See :mod:`lsdtt3._llm_status`.

    Raises
    ------
    ValueError
        If ``name`` is not a public export, or not a member of one.

    Authors
    -------
    S. M. Mudd
    """
    owner, _, member = name.partition(".")
    if owner not in namespace["__all__"]:
        raise ValueError(f"Unknown lsdtt3 public symbol: '{owner}'.")
    if member and not hasattr(namespace[owner], member):
        raise ValueError(f"'{owner}' has no member '{member}'.")
    if name in EXPERIMENTAL_OVERRIDES or owner in EXPERIMENTAL_OVERRIDES:
        return "experimental"
    return "ready" if owner in READY else "experimental"


def _without_authors(doc: str) -> str:
    """Drop the trailing ``Authors`` section, which is noise for a model."""
    return re.sub(r"\n*Authors\n-{3,}\n(?:[^\n]+(?:\n|$))+\s*$", "", doc).rstrip()


_EXPERIMENTAL_NOTE = ("[experimental] not in the guide's topics yet: it works, but its"
                      " interface may change and the guide does not teach it")


def _symbol_card(name: str, value: object, status: str = "ready") -> str:
    """Render one public symbol or member: its whole docstring, without the authors.

    A function's card starts with its typed signature, then its docstring
    (summary, parameters, return value); a class's docstring holds its
    constructor (or what returns it) and one line per property and method. A
    method (``Class.method``) or property card is that member's docstring. An
    experimental symbol's card says so on its second line.
    """
    doc = _without_authors(inspect.getdoc(value) or "signature unavailable")
    if not isinstance(value, (type, property)) and not re.match(r"^[A-Za-z_]\w*\(", doc):
        try:
            doc = f"{name.rpartition('.')[2]}{inspect.signature(value)}\n\n{doc}"
        except (TypeError, ValueError):
            pass
    note = f"{_EXPERIMENTAL_NOTE}\n" if status == "experimental" else ""
    return f"@{name}\n{note}{doc}\n"


def _resolve_symbol(name: str, namespace: Mapping[str, object]) -> object:
    """The object named ``name`` or ``Class.member``; ValueError if there is none."""
    owner, _, member = name.partition(".")
    if owner not in namespace["__all__"]:
        raise ValueError(
            f"Unknown lsdtt3 public symbol: '{name}'. "
            "Try helper_for_llm('search:<words>') to find it."
        )
    if not member:
        return namespace[owner]
    cls = namespace[owner]
    if member.startswith("_") or not hasattr(cls, member):
        raise ValueError(
            f"'{owner}' has no public member '{member}'. "
            f"Try helper_for_llm('symbol:{owner}') for its members."
        )
    value = inspect.getattr_static(cls, member)
    return value if isinstance(value, property) else getattr(cls, member)


def _words(text: str) -> list[str]:
    """Lower-case words of ``text``, splitting CamelCase and snake_case."""
    spaced = re.sub(r"([a-z0-9])([A-Z])", r"\1 \2", text)
    return re.findall(r"[a-z0-9]+", spaced.lower())


def _stem(word: str) -> str:
    """Fold a plain English plural, so "basins" finds "basin"."""
    return word[:-1] if len(word) > 3 and word.endswith("s") else word


def _first_line(value: object) -> str:
    doc = inspect.getdoc(value) or ""
    return doc.strip().splitlines()[0].strip() if doc.strip() else ""


def _search_entries(namespace: Mapping[str, object]):
    """Every public export and every public member of a class, as
    (name, name words, summary, docstring)."""
    for name in namespace["__all__"]:
        value = namespace[name]
        yield name, _words(name), _summary(value), inspect.getdoc(value) or ""
        if isinstance(value, type):
            for member, member_value, _is_property in _catalogue_members(value):
                yield (f"{name}.{member}", _words(member), _summary(member_value),
                       inspect.getdoc(member_value) or "")


def _search(query: str, namespace: Mapping[str, object], limit: int = 12) -> str:
    """Rank public exports and class members (``Class.method``) for a query.

    An entry matches if any query word is in its name or docstring; names are
    split at underscores and CamelCase, so a CLI parameter name such as
    ``run_low_concavity_detection`` matches the routine that documents it.
    Ranking, best first: the entry's own name is exactly the query words
    (``fill pits`` finds ``fill_pits``, ``basins`` finds ``FlowInfo.basins``);
    more distinct query words matched; more matched in the entry's own name
    (a method's name, not its class's); more matched in the first docstring
    line; fewer extra words in the name; then alphabetical. Experimental
    entries are marked.
    """
    words = {_stem(word) for word in _words(query)}
    scored = []
    for name, name_words, summary, doc in _search_entries(namespace):
        name_set = {_stem(word) for word in name_words}
        summary_set = {_stem(word) for word in _words(summary)}
        doc_set = {_stem(word) for word in _words(doc)} | summary_set | name_set
        matched = words & doc_set
        if not matched:
            continue
        exact = name_set == words
        scored.append((-exact, -len(matched), -len(words & name_set),
                       -len(words & summary_set), len(name_set - words), name.lower(),
                       name, summary))
    if not scored:
        return f"@search {query.strip()}\nnone\n"
    lines = [f"@search {query.strip()} ({len(scored)} matches, best {min(limit, len(scored))} shown;"
             " helper_for_llm('symbol:<name>') for one)"]
    for *_key, name, summary in sorted(scored)[:limit]:
        marker = "[experimental] " if _guide_status(name, namespace) == "experimental" else ""
        lines.append(f"{marker}{name}: {summary}")
    return "\n".join(lines) + "\n"


def _split_top_level(text: str) -> list[str]:
    """Split at commas outside brackets."""
    parts, depth, start = [], 0, 0
    for index, character in enumerate(text):
        if character in "([{":
            depth += 1
        elif character in ")]}":
            depth -= 1
        elif character == "," and depth == 0:
            parts.append(text[start:index])
            start = index + 1
    parts.append(text[start:])
    return [part.strip() for part in parts if part.strip()]


def _parameter_list(definition: str) -> str:
    """``a,b?,*,c?`` from ``def name(a: T, b: T = \\0, *, c=1) -> R``: names only,
    ``?`` marking a parameter with a default."""
    start = definition.find("(")
    if start < 0:
        return ""
    depth = 0
    for end in range(start, len(definition)):
        depth += definition[end] in "([{"
        depth -= definition[end] in ")]}"
        if depth == 0:
            break
    names = []
    for parameter in _split_top_level(definition[start + 1:end]):
        name = re.split(r"[:=]", parameter, maxsplit=1)[0].strip()
        if name in ("self", "/"):
            continue
        has_default = "=" in re.sub(r"\[.*?\]", "", parameter)
        names.append(name + ("?" if has_default and not name.startswith("*") else ""))
    return ",".join(names)


def _definitions(value: object) -> list[str]:
    """The signature texts of a callable: nanobind's, else inspect's."""
    native = inspect.unwrap(value)
    signatures = getattr(native, "__nb_signature__", None)
    if signatures:
        return [definition for definition, *_rest in signatures]
    try:
        return [f"def f{inspect.signature(value)}"]
    except (TypeError, ValueError):
        return []


def _summary(value: object) -> str:
    """The first docstring line that is not a signature."""
    for line in (inspect.getdoc(value) or "").splitlines():
        line = line.strip()
        if line and not re.match(r"^[A-Za-z_][\w.]*\(.*\)(\s*->.*)?$", line):
            return line
    return ""


def _catalogue_members(cls: type) -> list[tuple[str, object, bool]]:
    """Public members of a class as (name, value, is_property), as the API
    manifest lists them (tools/generate_public_api_manifest.py)."""
    members = []
    for name in sorted(name for name in dir(cls) if not name.startswith("_")):
        try:
            value = inspect.getattr_static(cls, name)
        except AttributeError:
            continue
        if isinstance(value, property):
            members.append((name, value, True))
            continue
        bound = getattr(cls, name, value)
        if _definitions(bound):
            members.append((name, bound, False))
    return members


def _catalogue_lines(name: str, value: object, namespace: Mapping[str, object]) -> list[str]:
    """One line for an export; for a class, one more per ready public member."""
    if not isinstance(value, type):
        definitions = _definitions(value)
        parameters = _parameter_list(definitions[0]) if definitions else "..."
        return [f"{name}({parameters}): {_summary(value)}"]
    init = _definitions(getattr(value, "__init__", None)) if "__init__" in vars(value) or \
        hasattr(inspect.unwrap(value.__init__), "__nb_signature__") else []
    parameters = _parameter_list(init[0]) if init else ""
    # no public constructor (made by another call): show the bare class name
    head = f"{name}({parameters})" if init and parameters != "*args,**kwargs" else name
    lines = [f"{head}: {_summary(value)}"]
    for member, member_value, is_property in _catalogue_members(value):
        if _guide_status(f"{name}.{member}", namespace) != "ready":
            continue
        if is_property:
            lines.append(f"  .{member}: {_summary(member_value)}")
        else:
            definitions = _definitions(member_value)
            lines.append(f"  .{member}({_parameter_list(definitions[0])}): {_summary(member_value)}")
    return lines


def _catalogue(topic: str, namespace: Mapping[str, object]) -> str:
    """The generated catalogue of a topic: every READY export whose topic it is."""
    names = _ready_names(topic, namespace)
    lines = [f"@catalogue {topic}: every ready export (name(args): summary; arg?=has a default;"
             f" helper_for_llm('symbol:<name>') for types and defaults)"]
    for name in names:
        lines.extend(_catalogue_lines(name, namespace[name], namespace))
    return "\n".join(lines) + "\n"


def _ready_names(topic: str, namespace: Mapping[str, object]) -> list[str]:
    return sorted((name for name in namespace["__all__"]
                   if guide_topic(name) == topic and _guide_status(name, namespace) == "ready"),
                  key=str.lower)


def _all(namespace: Mapping[str, object]) -> str:
    """The pasteable guide: the index, every topic's guidance, the names of all
    ready exports by topic, and the method names of the key classes."""
    parts = [_INDEX] + [_GUIDANCE[name] for name in _ALL_ORDER]
    names = ["@names every ready export, by topic (arguments: helper_for_llm('symbol:<name>'))"]
    names += [f"{topic}: {' '.join(_ready_names(topic, namespace))}" for topic in GUIDE_TOPICS]
    methods = ["@methods of the key classes (arguments: helper_for_llm('symbol:<Class>.<method>'))"]
    for cls in _ALL_KEY_CLASSES:
        members = [member for member, _value, _is_property in _catalogue_members(namespace[cls])
                   if _guide_status(f"{cls}.{member}", namespace) == "ready"]
        methods.append(f"{cls}: {' '.join(members)}")
    closing = ("@end Without tools, ask the user to run print(lsdtt3.helper_for_llm('symbol:<name>'))"
               " for the exact arguments of a name above, and paste the result.\n")
    return "".join(parts) + "\n".join(names) + "\n" + "\n".join(methods) + "\n" + closing


def _topic_card(topic: str, namespace: Mapping[str, object]) -> str:
    """A topic: its hand-written guidance, then its catalogue if it has one."""
    if topic in GUIDE_TOPICS:
        return _GUIDANCE[topic] + _catalogue(topic, namespace)
    return _GUIDANCE[topic]


def helper_for_llm(topic: str = "index", *, _namespace: Mapping[str, object]) -> str:
    """Return a compact API card selected by topic, symbol, or search query.

    Parameters
    ----------
    topic : str, default="index"
        Topic name, ``"symbol:<name>"``, ``"search:<words>"``, or ``"all"``.
    _namespace : collections.abc.Mapping[str, object]
        Public package namespace supplied by :mod:`lsdtt3`.

    Returns
    -------
    str
        Compact retrieval-oriented API context. It is a ``str`` subclass that
        notebooks display as plain text.

    Notes
    -----
    ``_namespace`` is supplied only by :mod:`lsdtt3` to avoid import cycles.
    It is intentionally private rather than part of the public helper API.

    Authors
    -------
    B. Gailleton
    S. M. Mudd
    """
    return _PlainText(_render(topic, _namespace))


def _render(topic: str, namespace: Mapping[str, object]) -> str:
    request = topic.strip().lower()
    if request in ("", "index", "help"):
        return _INDEX
    if request == "all":
        return _all(namespace)
    if request in ("catalogue", "catalog"):
        return _INDEX.splitlines()[0] + "\n" + "".join(
            _catalogue(name, namespace) for name in GUIDE_TOPICS)
    if request.startswith("symbol:"):
        name = topic.split(":", 1)[1].strip()
        value = _resolve_symbol(name, namespace)
        return _symbol_card(name, value, _guide_status(name, namespace))
    if request.startswith("search:"):
        return _search(topic.split(":", 1)[1], namespace)
    request = _ALIASES.get(request, request)
    if request not in _GUIDANCE:
        choices = ", ".join(sorted(_GUIDANCE))
        raise ValueError(f"Unknown lsdtt3 LLM topic: '{topic}'. Try: {choices}.")
    return _topic_card(request, namespace)
