"""Compact, retrieval-oriented API context for language models.

This module intentionally favours exact signatures and behavioural constraints
over prose.  It is not human-facing reference documentation.

Authors
-------
B. Gailleton
"""

from __future__ import annotations

from collections.abc import Mapping
import re


_INDEX = """LSDTT3/API1
T:R=Raster M=RasterMeta F=FlowInfo G=FlowGraph S=ChannelScope N=ChannelNetwork B=BasinCollection C=ChiCollection V=GeoDataFrame D=DataFrame
!: raster funcs pure; A[row,col]; nodata=R.mask_numpy(); V CRS=R projected metres; metres/years unless named otherwise.
@:terrain flow channel basin chi vector fetch hydraulic landscape cosmo drainage synthetic stats
?: helper_for_llm(T); symbol:<exact>; search:<words>; all
"""

_TOPICS = {
    "terrain": """@terrain
R.read(path)->R; R.from_numpy(a,mask,M)->R; R.write(path,format='cog'); R.to_numpy(); R.mask_numpy()
M.from_shape(rows,cols,x_min=0,y_min=0,cell_size=1,...)->M
resample_raster(R,target_meta,method='linear')->R; double_raster|halve_raster(R,method)->R
pad_to_nearest_power_of_two(R)->R; interpolate_to_nearest_power_of_two(R,direction='nearest',method='linear')->R
H=RasterHierarchy(R,storage='retain_allocated',keep_finest=True,levels=0); H.restrict_to(level,method='conservative_mean')->R; H.prolong_to(level,method='linear')->R
fill_pits(R,...)->R; breach_depressions(R,...)->R; slope(R)->R; multidirectional_hillshade(R,...)->R
polyfit_metrics(R,[metric],lambda_metres=0,...)->{metric:R}; classify_surface_shape(R|kmax,kmin,...)->{...}
spectral_analysis(R,...)->SpectralResult; frangi_filter(R,...)->FrangiResult; swath_profile(R,baseline:V,...)->SwathProfile
""",
    "flow": """@flow
F=FlowInfo(R,boundary_conditions='oooo',forced_outlets?,allow_pits_as_outlets?); F is immutable D8 topology.
F.drainage_area()->R; F.flow_direction()->R; F.flow_distance_to_outlet()->R; F.channel_scope(threshold_pixels=...)->S
F.scope_from_heads(heads)->S; F.channel_network(S)->N; F.basins(outlets?)->B; F.height_above_channel(R,S|N,...)->R
G=FlowGraph(R,method='mfd'|'dinf',boundary_conditions='oooo',...); G.drainage_area()->R
topographic_wetness_index(G,slope:R,...)->R; sediment_connectivity_index(...)->ConnectivityResult; soil_topographic_index(...)->SoilTopographicIndexResult
""",
    "channel": """@channel
extract_channel_heads(F|R,method=..., ...)->ChannelHeads; S=F.channel_scope(threshold_pixels=...)|F.scope_from_heads(heads)
N=F.channel_network(S); N.points_geodataframe()->V; N.lines_geodataframe()->V; N.links()->D; N.junctions()->D; N.junction_angles()->D
N.chi(R,m_over_n=.45,area_threshold?)->C; N.chi_sweep(R,m_over_n_values,...)->D
find_channel_candidates(R,...)->ChannelCandidatesResult; water_centrelines(...)->WaterCentrelineResult; repair_water_lines(...)->WaterLineRepairResult
""",
    "basin": """@basin
B=F.basins(outlets?); junction_basins(N,selection='all',outlets?)->JunctionBasinCollection
B.labels()->R; B.mask(id)->R; B.membership()->D; B.pixels()->D; B.polygons()->V; B.statistics(R,...)->D
basin_surface_metrics(B,R,...)->D; basin_channel_lengths(B,N,...)->D; basin_ridgeline_extrema(B,N,R,...)->V
divide_migration(R,F,N,...)->DivideMigrationResult; B.outline_ridgelines(...)->list[Ridgeline]
""",
    "chi": """@chi
C=N.chi(R,m_over_n=.45,area_threshold?); C.nodes()->D; C.raster()->R; C.fit(...)->ChiFit
ksn: use N.chi(...).fit/segment workflow; ChiInversion via N.chi(...), inversion helpers, or chi_ksn_inversion API.
N.chi_sweep(R,m_over_n_values,...)->D; chi_from_outlets/chi_from_outlets_inherited are F methods for outlet-defined chi.
! chi requires a D8 N/F and elevation R on the same grid; inspect N.points_geodataframe() for mapped network output.
""",
    "vector": """@vector
V inputs/outputs=GeoDataFrame; CRS must equal R projected metre CRS.
bounds_from_vector(V)->VectorBounds; rasterize_polygons(R,V,...)->R; clip_to_polygons(R,V)->R; sample_points(R,V)->V
burn_streams(R,V,...)->StreamBurnResult; water_centrelines(...)->WaterCentrelineResult; repair_water_lines(...)->WaterLineRepairResult
extract_valley_centrelines(R,N|S,...)->ValleyCentrelineResult; extract_valley_cross_sections(R,baseline:V,...)->ValleyCrossSections
""",
    "fetch": """@fetch
fetch_dem(bounds|V,dataset='cop30',...)->R; fetch_dem_for_outlets(outlets:V,...)->OutletFetchResult
bounds_from_vector(V)->VectorBounds; utm_epsg_for_point(x,y)->int; cop30_tiles_for_bbox(...)->list[str]
! fetch expects geographic bounds/outlets where documented; returned R carries projected grid metadata. Save explicitly: R.write(path).
""",
    "hydraulic": """@hydraulic
fill_hydraulic_surface(z:R,h:R,node_state:R,boundary_conditions='oooo',...)->HydraulicSurfaceFill
run_graphflood(z:R,node_state:R,rainfall_rate:R,friction_coefficient:R,...)->HydraulicResult
run_graphflood_optimised(z:R,node_state:R,rainfall_rate:R,friction_coefficient:R,...)->GraphFloodOptimisedResult
run_graphflood_optimised_v2(z:R,node_state:R,rainfall_rate:R,friction_coefficient:R,prescribed_qin:R|None=None,...)->GraphFloodOptimisedResult
run_graphflood_multigrid(z:R,node_state:R,rainfall_rate:R,friction_coefficient:R,backend='optimised_v2',...)->GraphFloodMultigridResult
run_graphflood_induced(z:R,node_state:R,rainfall_rate:R,friction_coefficient:R,qin_threshold,...)->GraphFloodInducedResult
node_state: 0=outside/nodata,1=active,2=ghost,3=outlet; ghosts retain initial depth and route prescribed_qin, only 3 drains; each active component needs an outlet unless check_valid_outlet=False.
! hydraulic routing is on z+h; stationary workflow fills h, not bed elevation.
""",
    "landscape": """@landscape
LandscapeModel(R,fluvial_method='d8'|'dinf',hillslope_model='linear'|'nonlinear'|'cellular',coupling_regime='auto'|..., ...)->L
L.step(dt_years)->LandscapeStepResult; L.run(n_steps,dt_years)->R; L.set_forcing(uplift_rate,baselevel_fall_rate=0)
hillslope_transport_explicit|implicit|adi|cellular(R,channels,...)->HillslopeTransportResult
route_sediment(R,F|G,...)->SedimentRoutingResult; trigger_landslides(R,...)->LandslideResult
""",
    "cosmo": """@cosmo
CRNColumn(): history owner; add_erosion_rate(rate_m_yr,duration_yr); add_burial(...)
crn_uniform_raster(R,value=1)->R; crn_production_scaling(R,...)->R; crn_steady_state(...); crn_step_change(...)
crn_transient_history(column,production_scaling:R,topographic_shielding:R?)->{'be10':R,'al26':R,'c14':R}
crn_invert_basin_erosion(...)->D; crn_outlet_concentration(...); crn_pair_diagnostics(...)
""",
    "drainage": """@drainage
find_internal_drainage(R,...)->InternalDrainageResult; depression_hierarchy(R,boundary_conditions='oooo',...)->DepressionHierarchyResult
result.labels->R; result.flow_direction->R; result.nodes->D (id,parent,spill,volume,...)
fill_pits_boundary_aware|fill_pits_with_outlets|fill_pits_with_sinks(R,...)->R
""",
    "synthetic": """@synthetic
synthetic_noise(M,...)->R; synthetic_mountain_range(M,...)->R; synthetic_hillslope_strip(M,...)->R
synthetic_network(M,...)->SyntheticNetworkResult; add_uniform_noise(R,...)->R
""",
    "stats": """@stats
linear_stats(...)->D; compare_bearings(...)->D; relief_ratio(R,...)->R
add|subtract|multiply|divide(R,R|scalar)->R; threshold_above|threshold_below(R,value)->R; mask_by_binary(R,mask)->R
""",
}

_ALIASES = {
    "raster": "terrain", "surface": "terrain", "channels": "channel",
    "networks": "channel", "basins": "basin", "vectors": "vector",
    "graphflood": "hydraulic", "hydraulics": "hydraulic", "lem": "landscape",
    "landslides": "landscape", "crn": "cosmo", "cosmogenic": "cosmo",
    "depression": "drainage", "internal": "drainage", "noise": "synthetic",
}

_CLASS_SIGNATURES = {
    "Raster": "Raster.read(path)->R | Raster.from_numpy(data,mask,meta)->R",
    "RasterMeta": "RasterMeta.from_shape(rows,cols,x_min=0,y_min=0,cell_size=1,...)->M",
    "FlowInfo": "FlowInfo(dem,boundary_conditions='oooo',forced_outlets?,allow_pits_as_outlets?)->F",
    "FlowGraph": "FlowGraph(dem,method='mfd'|'dinf',boundary_conditions='oooo',...)->G",
    "ChannelScope": "ChannelScope: F.channel_scope(...)|F.scope_from_heads(...) -> S",
    "ChannelNetwork": "ChannelNetwork: F.channel_network(S)->N",
    "BasinCollection": "BasinCollection: F.basins(outlets?)->B",
    "ChiCollection": "ChiCollection: N.chi(dem,m_over_n=.45,...)->C",
    "LandscapeModel": "LandscapeModel(dem,fluvial_method=...,hillslope_model=...,coupling_regime=...)->L",
    "CRNColumn": "CRNColumn(); add_erosion_rate(rate_m_yr,duration_yr); add_burial(...)",
}


def _symbol_card(name: str, value: object) -> str:
    """Render one symbol using only its first signature and summary lines."""
    target = getattr(value, "__wrapped__", value)
    doc = (getattr(target, "__doc__", None) or "").strip()
    lines = [line.strip() for line in doc.splitlines() if line.strip()]
    first = _CLASS_SIGNATURES.get(name, lines[0] if lines else "signature unavailable")
    summary = next((line for line in lines[1:] if not line.startswith("``")), "")
    return f"@{name}\n{first}" + (f"\n{summary}" if summary else "") + "\n"


def _search(query: str, namespace: Mapping[str, object]) -> str:
    """Return exact public symbols matching all case-insensitive query words."""
    words = set(re.findall(r"[a-z0-9]+", query.lower()))

    def symbol_words(name: str) -> set[str]:
        spaced = re.sub(r"([a-z0-9])([A-Z])", r"\1 \2", name)
        return set(re.findall(r"[a-z0-9]+", spaced.lower()))

    names = sorted(
        name for name in namespace["__all__"]
        if words <= symbol_words(name)
    )
    return "@search " + (" ".join(names) if names else "none") + "\n"


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
        Compact retrieval-oriented API context.

    Notes
    -----
    ``_namespace`` is supplied only by :mod:`lsdtt3` to avoid import cycles.
    It is intentionally private rather than part of the public helper API.

    Authors
    -------
    B. Gailleton
    """
    request = topic.strip().lower()
    if request in ("", "index", "help"):
        return _INDEX
    if request == "all":
        return _INDEX + "".join(_TOPICS[name] for name in sorted(_TOPICS))
    if request.startswith("symbol:"):
        name = topic.split(":", 1)[1].strip()
        if name not in _namespace["__all__"]:
            raise ValueError(f"Unknown lsdtt3 public symbol: '{name}'.")
        return _symbol_card(name, _namespace[name])
    if request.startswith("search:"):
        return _search(topic.split(":", 1)[1], _namespace)
    request = _ALIASES.get(request, request)
    if request not in _TOPICS:
        choices = ", ".join(sorted(_TOPICS))
        raise ValueError(f"Unknown lsdtt3 LLM topic: '{topic}'. Try: {choices}.")
    return _TOPICS[request]
