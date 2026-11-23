"""lsdtt-slope-stability's routines in Python: sinmap, FlowGraph.dinf_contour_width,
FlowInfo.rupture_planes -> RupturePlanes (failure_probability,
landslide_catalogues, slip_surfaces -> SlipSurfaces), landslide_frequency_area
and GroundwaterDupuit. The exact cases mirror the C++ drivers
(tests/cpp/test_slope_stability_*.cpp, test_groundwater_dupuit.cpp)."""

import math

import numpy as np
import pytest
from pyproj import CRS

import lsdtt3

pytestmark = [pytest.mark.unit]

WKT = CRS.from_epsg(32633).to_wkt()


def _raster(values, cell_size=10.0):
    values = np.asarray(values, dtype=np.float32)
    meta = lsdtt3.RasterMeta.from_shape(values.shape[0], values.shape[1], x_min=0.0, y_min=0.0,
                                        cell_size=cell_size, crs_wkt=WKT, epsg_code=32633)
    return lsdtt3.Raster.from_numpy(values, np.ones(values.shape, dtype=bool), meta)


def _east_plane(nr=12, nc=40, cell_size=10.0, gradient=0.5):
    """Elevation falls to the east; the last column is the channel."""
    z = 100.0 + gradient * cell_size * (nc - np.arange(nc))[None, :] * np.ones((nr, 1))
    return _raster(z, cell_size)


def _east_channel(nr=12, nc=40, cell_size=10.0):
    c = np.zeros((nr, nc))
    c[:, -1] = 1.0
    return _raster(c, cell_size)


def _hummocky(nr=40, nc=50, cell_size=10.0):
    r, c = np.mgrid[0:nr, 0:nc].astype(float)
    z = 400.0 - 4.0 * c + 15.0 * np.sin(c / 5.0) * np.cos(r / 6.0) + 8.0 * np.sin(r / 3.0 + c / 7.0)
    return _raster(z, cell_size)


def _valid(raster):
    return raster.mask_numpy().astype(bool)


# ------------------------------------------------------------------- SINMAP
@pytest.mark.exact_solution
def test_sinmap_factor_of_safety_bounds_on_a_plane():
    dem = _east_plane(gradient=0.4)
    graph = lsdtt3.FlowGraph(dem, method="dinf")
    out = lsdtt3.sinmap(graph.slope(), graph.drainage_area(), contour_width=graph.dinf_contour_width())
    S = 0.4
    theta = math.atan(S)
    a = out["specific_catchment_area"].to_numpy()
    rows, cols = slice(3, 9), slice(5, 30)
    a_inner = a[rows, cols].astype(float)
    # FS = (C + cos(theta) (1 - w r) t) / sin(theta), w = min(1, R a / (T sin)) (Pack et al. 1998).
    w_min = np.minimum(1.0, a_inner / (3000.0 * math.sin(theta)))   # (T/R)_max for FS_max
    w_max = np.minimum(1.0, a_inner / (2000.0 * math.sin(theta)))   # (T/R)_min for FS_min
    fs_min = (0.0 + math.cos(theta) * (1.0 - 0.5 * w_max) * math.tan(math.radians(30))) / math.sin(theta)
    fs_max = (0.25 + math.cos(theta) * (1.0 - 0.5 * w_min) * 1.0) / math.sin(theta)
    np.testing.assert_allclose(out["fs_min"].to_numpy()[rows, cols], fs_min, rtol=1e-5)
    np.testing.assert_allclose(out["fs_max"].to_numpy()[rows, cols], fs_max, rtol=1e-5)
    table = out["class_areas"]
    assert list(table.columns) == ["region_id", "stability_class", "n_cells", "area_m2", "fraction_of_region"]
    assert table["fraction_of_region"].sum() == pytest.approx(1.0)
    assert out["n_evaluated"] == int(table["n_cells"].sum())


@pytest.mark.exact_solution
def test_dinf_contour_width_axis_and_diagonal():
    east = lsdtt3.FlowGraph(_east_plane(), method="dinf").dinf_contour_width().to_numpy()
    assert east[5, 10] == pytest.approx(10.0)
    r, c = np.mgrid[0:20, 0:20].astype(float)
    diag = lsdtt3.FlowGraph(_raster(500.0 - 3.0 * (r + c)), method="dinf").dinf_contour_width().to_numpy()
    assert diag[8, 8] == pytest.approx(10.0 / math.sqrt(2.0), rel=1e-5)
    with pytest.raises(ValueError, match="dinf"):
        lsdtt3.FlowGraph(_east_plane(), method="mfd").dinf_contour_width()


@pytest.mark.property
def test_sinmap_slope_area_curves_names_and_saturation_line():
    df = lsdtt3.sinmap_slope_area_curves(1.0, 10.0, 1.0e4)
    assert list(df.columns) == ["curve", "slope", "specific_catchment_area_m"]
    assert {"always_saturated_above", "never_saturated_below"} <= set(df["curve"])
    line = df[df["curve"] == "never_saturated_below"]
    # a = sin(theta) / x2, x2 = 1 / (T/R)_max ... the CLI's (T/R)_min -> x2 coding: a = (T/R) sin(theta).
    s = line["slope"].to_numpy()
    a = line["specific_catchment_area_m"].to_numpy()
    sin = s / np.sqrt(1.0 + s * s)
    ratio = a / sin
    assert np.allclose(ratio, ratio[0], rtol=1e-9) and ratio[0] in (pytest.approx(2000.0), pytest.approx(3000.0))


@pytest.mark.metamorphic
@pytest.mark.negative
def test_sinmap_regions_and_refusals():
    dem = _east_plane(gradient=0.4)
    graph = lsdtt3.FlowGraph(dem, method="dinf")
    reg = np.zeros((12, 40))
    reg[6:, :] = 1
    region = _raster(reg)
    out = lsdtt3.sinmap(graph.slope(), graph.drainage_area(), region=region,
                        region_parameters={0: {}, 1: {"sinmap_dimensionless_cohesion_min": 0.2,
                                                      "rupture_cohesion_pa": 5.0}})
    fs = out["fs_min"].to_numpy()
    assert fs[8, 10] > fs[3, 10]
    with pytest.raises(ValueError, match="region_parameters"):
        lsdtt3.sinmap(graph.slope(), graph.drainage_area(), region=region)
    with pytest.raises(ValueError, match="unknown parameter 'sinmap_cohesion'"):
        lsdtt3.sinmap(graph.slope(), graph.drainage_area(), region=region,
                      region_parameters={0: {"sinmap_cohesion": 1.0}, 1: {}})


# ------------------------------------------------------------ rupture planes
@pytest.mark.property
@pytest.mark.metamorphic
def test_rupture_planes_on_a_plane_and_threads():
    dem = _east_plane()
    flow = lsdtt3.FlowInfo(dem)
    rp = flow.rupture_planes(dem, channel=_east_channel())
    sf = rp.safety_factor.to_numpy()
    ok = _valid(rp.safety_factor)
    assert not ok[:, -1].any()                      # channel cells are not sources
    # The column next to the channel has no daylight candidate below it.
    assert ok[:, :-2].all() and not ok[:, -2].any()
    assert rp.n_channel == 12 and rp.n_sources == 12 * 38
    unstable = rp.unstable.to_numpy()
    assert np.array_equal(unstable[ok] == 1, sf[ok] <= 1.0)
    assert rp.n_unstable == int((unstable[ok] == 1).sum())
    assert rp.daylight.dtype == np.int64 and rp.daylight.shape == (12 * 40,)
    rp8 = flow.rupture_planes(dem, channel=_east_channel(), n_threads=8)
    assert np.array_equal(rp8.safety_factor.to_numpy()[ok], sf[ok])
    # Higher cohesion is never less stable.
    rp_c = flow.rupture_planes(dem, channel=_east_channel(), cohesion_pa=5000.0)
    assert (rp_c.safety_factor.to_numpy()[ok] >= sf[ok]).all()


@pytest.mark.metamorphic
def test_rupture_channel_from_threshold_matches_contributing_pixels():
    dem = _hummocky()
    flow = lsdtt3.FlowInfo(dem)
    rp = flow.rupture_planes(dem, threshold_contributing_pixels=50)
    npix = flow.contributing_pixels().to_numpy()
    assert np.array_equal(rp.channel.to_numpy() == 1, npix >= 50)


@pytest.mark.exact_solution
def test_zero_width_uncertainty_reproduces_the_deterministic_run():
    dem = _hummocky()
    flow = lsdtt3.FlowInfo(dem)
    rp = flow.rupture_planes(dem, threshold_contributing_pixels=50, cohesion_pa=800.0, friction_coefficient=0.6)
    pf = rp.failure_probability(cohesion_pa_min=800.0, cohesion_pa_max=800.0, friction_coefficient_min=0.6,
                                friction_coefficient_max=0.6, uncertainty_samples_per_parameter=3)
    ok = _valid(rp.safety_factor)
    assert pf["n_samples"] == 1
    np.testing.assert_array_equal(pf["failure_probability"].to_numpy()[ok], rp.unstable.to_numpy()[ok])
    wide = rp.failure_probability(uncertainty_samples_per_parameter=4)
    assert wide["n_samples"] == 16
    p = wide["failure_probability"].to_numpy()[ok]
    assert ((p >= 0.0) & (p <= 1.0)).all()


@pytest.mark.property
def test_landslide_catalogues_match_their_labels():
    dem = _hummocky()
    flow = lsdtt3.FlowInfo(dem)
    rp = flow.rupture_planes(dem, threshold_contributing_pixels=50, cohesion_pa=200.0, friction_coefficient=0.4)
    assert rp.n_unstable > 0
    cat = rp.landslide_catalogues(min_daylight_neighbours=0)
    by = cat["by_daylight"]
    assert list(by.columns[:4]) == ["landslide_id", "daylight_x", "daylight_y", "n_cells"]
    assert "n_pieces" in by.columns and "n_daylight_points" in cat["amalgamated"].columns
    lab = cat["label"].to_numpy()
    lab_ok = _valid(cat["label"])
    counts = {i: int(((lab == i) & lab_ok).sum()) for i in by["landslide_id"]}
    assert all(counts[i] == n for i, n in zip(by["landslide_id"], by["n_cells"]))
    np.testing.assert_allclose(by["area_m2"], by["n_cells"] * 100.0)
    assert cat["amalgamated"]["n_cells"].sum() == cat["n_landslide_cells"] == by["n_cells"].sum()


@pytest.mark.metamorphic
@pytest.mark.negative
def test_region_parameters_and_refusals():
    dem = _hummocky()
    flow = lsdtt3.FlowInfo(dem)
    reg = np.zeros((40, 50))
    reg[:, 25:] = 7
    region = _raster(reg)
    rp = flow.rupture_planes(dem, threshold_contributing_pixels=50, region=region,
                             region_parameters={0: {}, 7: {"rupture_cohesion_pa": 9000.0}})
    base = flow.rupture_planes(dem, threshold_contributing_pixels=50)
    ok = _valid(rp.safety_factor)
    right = np.zeros_like(ok)
    right[:, 25:] = True
    assert (rp.safety_factor.to_numpy()[ok & right] > base.safety_factor.to_numpy()[ok & right]).all()
    np.testing.assert_array_equal(rp.safety_factor.to_numpy()[ok & ~right], base.safety_factor.to_numpy()[ok & ~right])
    with pytest.raises(ValueError, match="region_parameters"):
        flow.rupture_planes(dem, region=region)
    with pytest.raises(ValueError, match="needs a region raster"):
        flow.rupture_planes(dem, region_parameters={0: {}})
    with pytest.raises(ValueError, match="normal"):
        base.failure_probability(cohesion_distribution="lognormal")


# ------------------------------------------------------------ slip surfaces
@pytest.mark.property
@pytest.mark.parametrize("family", ["plane_integrated", "circular", "bilinear", "structural", "ellipsoid"])
def test_slip_surface_families_run_and_report(family):
    dem = _hummocky()
    flow = lsdtt3.FlowInfo(dem)
    rp = flow.rupture_planes(dem, threshold_contributing_pixels=50, cohesion_pa=500.0, friction_coefficient=0.5)
    ss = rp.slip_surfaces(family=family, screen_safety_factor=3.0)
    assert ss.family == family
    assert ss.method == {"plane_integrated": "ordinary", "circular": "bishop", "bilinear": "janbu",
                         "structural": "janbu", "ellipsoid": "xie_revised"}[family]
    assert ss.n_evaluated + ss.n_screened <= rp.n_sources
    ev = ss.evaluated.to_numpy()
    assert int((ev[_valid(ss.evaluated)] == 1).sum()) == ss.n_evaluated
    sf = ss.safety_factor.to_numpy()
    ok = _valid(ss.safety_factor)
    assert ss.n_unstable == int((sf[ok] <= 1.0).sum())
    cats = ss.body_catalogues()
    assert {"by_toe", "amalgamated", "label_amalgamated"} <= set(cats)
    pf = ss.failure_probability(uncertainty_samples_per_parameter=2)
    assert pf["n_samples"] == 4
    if family == "ellipsoid":
        assert "cover_failure_probability" in pf
        with pytest.raises(ValueError, match="2-D"):
            ss.spencer_check()
    else:
        sp = ss.spencer_check()
        assert isinstance(sp["safety_factor"], lsdtt3.Raster) and sp["n_no_solution"] >= 0


@pytest.mark.negative
def test_slip_surface_refusals():
    dem = _hummocky()
    rp = lsdtt3.FlowInfo(dem).rupture_planes(dem, threshold_contributing_pixels=50)
    with pytest.raises(ValueError, match="family"):
        rp.slip_surfaces(family="wedge")
    with pytest.raises(ValueError, match="method"):
        rp.slip_surfaces(method="morgenstern_price")
    with pytest.raises(ValueError, match="ellipsoid family only"):
        rp.slip_surfaces(family="circular", layer_depths=[_raster(np.full((40, 50), 2.0))])


# ------------------------------------------------------------ frequency - area
@pytest.mark.exact_solution
def test_landslide_frequency_area_recovers_a_power_law_exponent():
    rng = np.random.default_rng(3)
    alpha = 2.4
    cells = np.floor(5.0 * (1.0 - rng.random(20000)) ** (-1.0 / (alpha - 1.0))).astype(float)
    density, fit = lsdtt3.landslide_frequency_area(cells * 100.0, 100.0, min_tail=200)
    assert list(density.columns) == ["area_lower_m2", "area_upper_m2", "area_centre_m2", "n_landslides",
                                     "probability_density_per_m2"]
    assert fit["method"] == "mle_discrete_ks_xmin" and fit["n_landslides"] == 20000
    assert fit["eta"] == pytest.approx(-alpha, abs=0.1)
    assert fit["n_tail"] >= 200
    with pytest.raises(ValueError, match="cell_area_m2"):
        lsdtt3.landslide_frequency_area([100.0], 0.0)


# ------------------------------------------------------------------ groundwater
@pytest.mark.exact_solution
def test_groundwater_kinematic_steady_state_on_a_steep_plane():
    K, f, S, d, cs = 1e-4, 2e-8, 0.5, 2.0, 10.0
    nr, nc = 5, 101
    dem = _east_plane(nr=nr, nc=nc, cell_size=cs, gradient=S)
    aq = lsdtt3.GroundwaterDupuit(dem, _east_channel(nr, nc, cs), regolith_thickness_metres=d,
                                  hydraulic_conductivity_m_per_s=K)
    rep = aq.run_to_steady(f, steady_tolerance=1e-9)
    assert rep["steady_method"] == "implicit" and rep["steps"] >= 1
    h = aq.saturated_thickness().to_numpy()[2]
    theta = math.atan(S)
    x = (np.arange(nc) + 1.0) * cs
    exact = f * x / (K * math.sin(theta) * math.cos(theta))
    np.testing.assert_allclose(h[20:nc - 10], exact[20:nc - 10], rtol=0.01)
    rates = aq.steady_rates()
    assert abs(rates["imbalance_m3_per_s"]) <= 1e-9 * rates["recharge_m3_per_s"]
    np.testing.assert_allclose(aq.water_table_depth().to_numpy()[2, :-1], d - h[:-1], atol=1e-5)


@pytest.mark.metamorphic
@pytest.mark.property
def test_groundwater_implicit_and_explicit_agree_and_budget_closes():
    dem = _hummocky(nr=20, nc=25)
    flow = lsdtt3.FlowInfo(dem)
    chan = flow.rupture_planes(dem, threshold_contributing_pixels=30).channel
    kw = dict(regolith_thickness_metres=2.0, hydraulic_conductivity_m_per_s=1e-4)
    a = lsdtt3.GroundwaterDupuit(dem, chan, **kw)
    b = lsdtt3.GroundwaterDupuit(dem, chan, **kw)
    a.run_to_steady(3e-8, steady_tolerance=1e-10)
    b.run_to_steady(3e-8, steady_method="explicit", steady_tolerance=1e-10)
    np.testing.assert_allclose(a.saturated_thickness().to_numpy(), b.saturated_thickness().to_numpy(),
                               atol=1e-6, equal_nan=True)
    # A storm: both steppers close the budget.
    for method in ("explicit", "implicit"):
        g = lsdtt3.GroundwaterDupuit(dem, chan, **kw)
        g.run_to_steady(3e-8)
        g.run(86400.0, 1e-6, transient_method=method, implicit_max_dt_s=3600.0)
        bud = g.budget()
        lhs = bud["recharge_m3"] + bud["clipped_added_m3"]
        rhs = bud["storage_change_m3"] + bud["outflow_fixed_head_m3"] + bud["seepage_m3"]
        assert lhs == pytest.approx(rhs, rel=1e-8), method
        assert bud["seepage_m3"] > 0.0


@pytest.mark.metamorphic
def test_a_water_table_lowers_the_rupture_safety_factor():
    dem = _hummocky()
    flow = lsdtt3.FlowInfo(dem)
    dry = flow.rupture_planes(dem, threshold_contributing_pixels=50)
    aq = lsdtt3.GroundwaterDupuit(dem, dry.channel, regolith_thickness_metres=3.0, hydraulic_conductivity_m_per_s=1e-5)
    aq.run_to_steady(1e-7)
    wet = flow.rupture_planes(dem, threshold_contributing_pixels=50, water_table_depth=aq.water_table_depth())
    ok = _valid(dry.safety_factor) & _valid(wet.safety_factor)
    assert (wet.safety_factor.to_numpy()[ok] <= dry.safety_factor.to_numpy()[ok] * (1 + 1e-6)).all()
    assert wet.n_unstable >= dry.n_unstable


@pytest.mark.negative
def test_groundwater_refusals():
    dem = _east_plane()
    chan = _east_channel()
    with pytest.raises(ValueError, match="transmissivity_profile"):
        lsdtt3.GroundwaterDupuit(dem, chan, transmissivity_profile="linear")
    aq = lsdtt3.GroundwaterDupuit(dem, chan)
    with pytest.raises(ValueError, match="steady_method"):
        aq.run_to_steady(1e-8, steady_method="newton")
    with pytest.raises(ValueError, match="transient_method"):
        aq.run(3600.0, 1e-8, transient_method="crank_nicolson")


# ------------------------------------------------------- DEM of difference
def _change_pair(n=160, cell=30.0, seed=3):
    rng = np.random.default_rng(seed)
    r, c = np.mgrid[0:n, 0:n].astype(float)
    dem = _raster(1000.0 + 15.0 * c + 3.0 * r, cell)
    dh = rng.normal(0.0, 1.0, (n, n))
    dh[40:60, 40:70] = -25.0
    dh[100:110, 100:120] = 20.0
    return dem, _raster(dh, cell), dh


def _hazard_profile(cell_size=10.0):
    """Twenty steep columns, five that carry debris, then a gentle reach."""
    gradient = np.r_[np.full(20, 1.0), np.full(5, 0.5), np.full(14, 0.05)]
    z = np.r_[np.cumsum((gradient * cell_size)[::-1])[::-1], 0.0]
    return _raster(np.tile(z, (20, 1)), cell_size)


@pytest.mark.exact_solution
def test_coseismic_hazard_area_exact_profile():
    out = lsdtt3.coseismic_hazard_area(_hazard_profile(), boundary_conditions="popo")
    a = out["hazard_area"].to_numpy()
    expected = np.r_[10.0 * np.arange(1, 21), np.full(6, 200.0), np.zeros(14)]
    np.testing.assert_allclose(a, np.tile(expected, (20, 1)), atol=1e-3)
    assert out["n_source_cells"] == 400
    assert out["n_zero_hazard_cells"] == 280
    assert out["source_area_m2"] == 40000.0


@pytest.mark.negative
def test_coseismic_hazard_area_refuses_inverted_thresholds():
    with pytest.raises(ValueError, match="stopping slope"):
        lsdtt3.coseismic_hazard_area(_hazard_profile(), initiation_slope=0.2, stopping_slope=0.5)


@pytest.mark.exact_solution
def test_dem_change_finds_a_planted_scar():
    dem, difference, _ = _change_pair()
    out = lsdtt3.dem_change_landslides(difference, dem, detection_sigma=4.5, seed_sigma=5.0, n_null_simulations=2,
                                       n_threads=1)
    patches = out["patches"]
    assert len(patches) == 1 and patches.change[0] == "loss" and patches.n_cells[0] == 600
    assert patches.area_m2[0] == pytest.approx(600 * 900.0)
    assert patches.mean_change_m[0] == pytest.approx(-25.0, abs=0.1)
    assert out["uncertainty"]["nmad_metres"] == pytest.approx(1.0, abs=0.08)
    assert out["change_class"].to_numpy()[50, 50] == -1
    assert out["false_detection"]["n_simulations"] == 2


@pytest.mark.property
def test_dem_change_sign_convention_and_stable_mask():
    dem, difference, dh = _change_pair()
    flipped = _raster(-dh, 30.0)
    a = lsdtt3.dem_change_landslides(difference, dem, detection_sigma=4.5, seed_sigma=5.0, n_null_simulations=0)
    b = lsdtt3.dem_change_landslides(flipped, dem, is_newer_minus_older=False, detection_sigma=4.5, seed_sigma=5.0,
                                     n_null_simulations=0)
    assert a["patches"].n_cells.tolist() == b["patches"].n_cells.tolist()
    stable = np.ones(dh.shape)
    stable[35:65, 35:75] = 0.0
    stable[95:115, 95:125] = 0.0
    c = lsdtt3.dem_change_landslides(difference, dem, stable_mask=_raster(stable, 30.0), detection_sigma=4.5,
                                     seed_sigma=5.0, n_null_simulations=0)
    assert c["uncertainty"]["used_stable_mask"] is True


@pytest.mark.negative
def test_dem_change_refuses_a_bad_sign():
    dem, difference, _ = _change_pair()
    with pytest.raises(ValueError):
        lsdtt3.dem_change_landslides(difference, dem, change_sign="up")
