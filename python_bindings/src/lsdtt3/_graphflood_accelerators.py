"""Higher-level CPU acceleration strategies for GraphFlood."""

from collections import deque
from dataclasses import dataclass
from typing import Any, Callable, Sequence

import numpy as np

from ._compat import MISSING, resolve_deprecated_keyword


def _normalize_convergence_options(
    absolute: float | None,
    relative: float | None,
    options: dict[str, Any],
) -> dict[str, Any]:
    options = dict(options)
    pairs = (
        ("median_dhw_convergence_tolerance_metres", "tol_median_abs_dhw", absolute),
        ("median_dhw_relative_convergence_tolerance", "tol_median_rel_dhw", relative),
    )
    for canonical, legacy, explicit in pairs:
        canonical_value = MISSING if explicit is None else explicit
        if canonical in options:
            if canonical_value is not MISSING:
                raise TypeError(f"{canonical!r} was supplied twice")
            canonical_value = options.pop(canonical)
        value = resolve_deprecated_keyword(
            canonical, canonical_value, legacy, options.pop(legacy, MISSING),
            stacklevel=3,
        )
        if value is not MISSING:
            options[canonical] = value
    return options


@dataclass(frozen=True)
class GraphFloodMultigridResult:
    """Final result plus every coarse-to-fine level.

    Parameters
    ----------
    result, level_results, factors
        Values returned by :func:`run_graphflood_multigrid`.

    Attributes
    ----------
    result : object
        Result from the full-resolution level.
    level_results : tuple
        Backend result at each requested factor.
    factors : tuple of int
        Coarsening factors in execution order.

    Authors
    -------
    B. Gailleton
    """

    result: Any
    level_results: tuple[Any, ...]
    factors: tuple[int, ...]

    def __getattr__(self, name: str) -> Any:
        return getattr(self.result, name)


@dataclass(frozen=True)
class GraphFloodInducedResult:
    """Reduced-domain result and prepass diagnostics.

    Parameters
    ----------
    result, prepass, active_mask, active_fraction, expansion_results
        Values returned by :func:`run_graphflood_induced`.

    Attributes
    ----------
    result : object
        Final reduced-domain backend result.
    prepass : object
        One-step full-domain routing result.
    active_mask : numpy.ndarray
        Final retained hydraulic domain.
    active_fraction : float
        Retained fraction of the original active cells.
    expansion_results : tuple
        Result from each domain-expansion solve.
    source_discharge, outlet_qin, relative_mass_error : float
        Full forcing, final outlet accumulation, and their relative mismatch.

    Authors
    -------
    B. Gailleton
    """

    result: Any
    prepass: Any
    active_mask: np.ndarray
    active_fraction: float
    expansion_results: tuple[Any, ...]
    source_discharge: float
    outlet_qin: float
    relative_mass_error: float

    def __getattr__(self, name: str) -> Any:
        return getattr(self.result, name)


def _backend(backend: str | Callable[..., Any]) -> tuple[Callable[..., Any], str]:
    if callable(backend):
        return backend, "custom"
    from . import run_graphflood, run_graphflood_optimised, run_graphflood_optimised_v2

    functions = {
        "reference": run_graphflood,
        "optimised": run_graphflood_optimised,
        "optimised_v2": run_graphflood_optimised_v2,
    }
    try:
        return functions[backend], backend
    except KeyError as error:
        raise ValueError(
            "backend must be 'reference', 'optimised', 'optimised_v2', or callable"
        ) from error


@dataclass(frozen=True)
class _MultigridTransferPolicy:
    """Field-specific transfer choices; kept explicit at the driver seam."""

    static_method: str = "conservative_mean"
    depth_restrict_method: str = "conservative_mean"
    depth_prolong_method: str = "linear"


_TRANSFER = _MultigridTransferPolicy()


def _same_grid(reference: Any, candidate: Any) -> bool:
    a, b = reference.metadata, candidate.metadata
    return (
        reference.shape == candidate.shape
        and a.x_min == b.x_min and a.y_min == b.y_min
        and a.cell_size == b.cell_size and a.crs_wkt == b.crs_wkt
    )


def _require_field(reference: Any, candidate: Any, name: str, active: np.ndarray) -> None:
    if not _same_grid(reference, candidate):
        raise ValueError(f"{name} must share the bed grid")
    valid = candidate.mask_numpy()
    if not np.all(valid[active]):
        raise ValueError(f"{name} must be valid on every active hydraulic cell")
    values = candidate.to_numpy()
    if not np.isfinite(values[active]).all():
        raise ValueError(f"{name} must be finite on every active hydraulic cell")


def _validate_multigrid_inputs(
    bed: Any, node_state: Any, rainfall_rate: Any,
    friction_coefficient: Any, initial_water_depth: Any | None,
) -> np.ndarray:
    if not _same_grid(bed, node_state):
        raise ValueError("node_state must share the bed grid")
    states = node_state.to_numpy()
    active = bed.mask_numpy() & node_state.mask_numpy() & (states != 0)
    _require_field(bed, rainfall_rate, "rainfall_rate", active)
    _require_field(bed, friction_coefficient, "friction_coefficient", active)
    if initial_water_depth is not None:
        _require_field(bed, initial_water_depth, "initial_water_depth", active)
    return active


def _coarse_node_state(bed: Any, node_state: Any, factor: int) -> Any:
    from . import Raster

    rows, cols = bed.shape
    states = node_state.to_numpy()
    active = bed.mask_numpy() & node_state.mask_numpy() & (states != 0)
    coarse_shape = (rows // factor, cols // factor)
    blocks = states.reshape(coarse_shape[0], factor, coarse_shape[1], factor)
    active_blocks = active.reshape(coarse_shape[0], factor, coarse_shape[1], factor)
    coarse_active = active_blocks.any(axis=(1, 3))
    coarse_outlet = ((blocks == 3) & active_blocks).any(axis=(1, 3))
    coarse_states = np.where(coarse_outlet, 3.0, np.where(coarse_active, 1.0, 0.0))
    meta = bed.metadata
    meta = type(meta).from_shape(
        *coarse_shape, x_min=meta.x_min, y_min=meta.y_min,
        cell_size=meta.cell_size * factor, nodata_value=meta.nodata_value,
        crs_wkt=meta.crs_wkt, epsg_code=meta.epsg_code,
    )
    return Raster.from_numpy(coarse_states.astype(np.float32), coarse_active, meta)


def _coarse_points(points: Any, factor: int, shape: tuple[int, int]) -> np.ndarray:
    if points is None:
        return np.empty((0, 3), dtype=np.float64)
    points = np.asarray(points, dtype=np.float64)
    if points.size == 0:
        return np.empty((0, 3), dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("point_inflows must be an N x 3 array")
    grouped: dict[tuple[int, int], float] = {}
    for row, col, discharge in points:
        if (
            not np.isfinite(row) or not np.isfinite(col) or not np.isfinite(discharge)
            or row != np.floor(row) or col != np.floor(col) or discharge < 0.0
        ):
            raise ValueError(
                "point_inflows must contain integral coordinates and non-negative discharge"
            )
        if row < 0 or row >= shape[0] or col < 0 or col >= shape[1]:
            raise ValueError("point_inflows contains a point outside the bed grid")
        key = (int(row) // factor, int(col) // factor)
        grouped[key] = grouped.get(key, 0.0) + float(discharge)
    return np.asarray([(r, c, q) for (r, c), q in grouped.items()], dtype=np.float64)


def run_graphflood_multigrid(
    bed: Any,
    node_state: Any,
    rainfall_rate: Any,
    friction_coefficient: Any,
    *,
    backend: str | Callable[..., Any] = "optimised_v2",
    factors: Sequence[int] = (8, 4, 2, 1),
    steps_per_level: int | Sequence[int] = 100,
    point_inflows: Any = None,
    initial_water_depth: Any = None,
    median_dhw_convergence_tolerance_metres: float | None = None,
    median_dhw_relative_convergence_tolerance: float | None = None,
    **backend_options: Any,
) -> GraphFloodMultigridResult:
    """Run a compatible GraphFlood backend from coarse to fine grids.

    Parameters
    ----------
    bed, node_state, rainfall_rate, friction_coefficient : Raster
        Full-resolution hydraulic inputs.
    backend : str or callable, default="optimised_v2"
        Bundled backend name or compatible GraphFlood function.
    factors : sequence of int, default=(8, 4, 2, 1)
        Strictly descending integer refinement factors ending at one.
    steps_per_level : int or sequence of int, default=100
        Solver iterations allocated to each level.
    point_inflows : array-like, optional
        Rows of ``(row, col, discharge_m3_s)``.
    initial_water_depth : Raster, optional
        Full-resolution initial depth.
    **backend_options
        Options forwarded to every backend call.

    Returns
    -------
    GraphFloodMultigridResult
        Final result and per-level results.

    Authors
    -------
    B. Gailleton
    """
    from . import Raster, RasterHierarchy

    backend_options = _normalize_convergence_options(
        median_dhw_convergence_tolerance_metres,
        median_dhw_relative_convergence_tolerance,
        backend_options,
    )

    solver, backend_name = _backend(backend)
    factors = tuple(int(value) for value in factors)
    if not factors or factors[-1] != 1 or any(value <= 0 for value in factors):
        raise ValueError("factors must be positive and end with 1")
    if any(a <= b or a % b != 0 for a, b in zip(factors, factors[1:])):
        raise ValueError("factors must be strictly descending integer multiples")
    if isinstance(steps_per_level, int):
        steps = (steps_per_level,) * len(factors)
    else:
        steps = tuple(int(value) for value in steps_per_level)
        if len(steps) != len(factors):
            raise ValueError("steps_per_level must match factors")

    if any(bed.shape[0] % factor or bed.shape[1] % factor for factor in factors):
        raise ValueError("each coarsening factor must divide both raster dimensions")

    active = _validate_multigrid_inputs(
        bed, node_state, rainfall_rate, friction_coefficient, initial_water_depth
    )
    # The built-in v2 path keeps all hierarchy transfers and solver hand-offs
    # in C++. Other backends deliberately retain the transparent Python driver.
    if backend_name == "optimised_v2":
        from ._core import _run_graphflood_optimised_v2_multigrid_native

        results, _ = _run_graphflood_optimised_v2_multigrid_native(
            bed, node_state, rainfall_rate, friction_coefficient,
            list(factors), list(steps), point_inflows=point_inflows,
            initial_water_depth=initial_water_depth, **backend_options,
        )
        results = tuple(results)
        return GraphFloodMultigridResult(results[-1], results, factors)

    zero_depth = Raster.from_numpy(
        np.zeros(bed.shape, dtype=np.float32), active, bed.metadata
    )
    depth_source = zero_depth if initial_water_depth is None else initial_water_depth
    if depth_source.shape != bed.shape:
        raise ValueError("initial_water_depth must match bed")

    # Static fields are restricted once from the finest data and retained.
    bed_h = RasterHierarchy(bed, storage="retain_allocated")
    friction_h = RasterHierarchy(friction_coefficient, storage="retain_allocated")
    depth_h = RasterHierarchy(depth_source, storage="retain_allocated")

    # Rainfall is averaged over the complete coarse cell. Invalid source cells
    # contribute zero rainfall but remain excluded by the coarse node state.
    rainfall_values = np.where(active, rainfall_rate.to_numpy(), 0.0).astype(np.float32)
    rainfall_transfer = Raster.from_numpy(
        rainfall_values, np.ones(bed.shape, dtype=bool), rainfall_rate.metadata
    )
    rainfall_h = RasterHierarchy(rainfall_transfer, storage="retain_allocated")

    factor_levels = {1: 0}
    for factor in sorted(set(factors)):
        if factor != 1:
            indices = [
                hierarchy.add_level_for_factor(factor)
                for hierarchy in (bed_h, friction_h, depth_h, rainfall_h)
            ]
            if len(set(indices)) != 1:
                raise RuntimeError("GraphFlood hierarchy layouts are inconsistent")
            factor_levels[factor] = indices[0]

    # Materialise static inputs and initial water depth by restriction only.
    materialise_factors = sorted((factor for factor in factors if factor != 1))
    for factor in materialise_factors:
        target = factor_levels[factor]
        for hierarchy, method in (
            (bed_h, _TRANSFER.static_method),
            (friction_h, _TRANSFER.static_method),
            (rainfall_h, _TRANSFER.static_method),
            (depth_h, _TRANSFER.depth_restrict_method),
        ):
            hierarchy.restrict_to(target, method)

    level_results = []
    for factor, max_steps in zip(factors, steps):
        level = factor_levels[factor]
        level_state = node_state if factor == 1 else _coarse_node_state(bed, node_state, factor)
        level_active = level_state.mask_numpy()
        source_bed = bed_h.level(level)
        source_friction = friction_h.level(level)
        source_rain = rainfall_h.level(level)
        level_bed = Raster.from_numpy(source_bed.to_numpy(), level_active, source_bed.metadata)
        level_friction = Raster.from_numpy(
            source_friction.to_numpy(), level_active, source_friction.metadata
        )
        level_rain = Raster.from_numpy(source_rain.to_numpy(), level_active, source_rain.metadata)

        if level != depth_h.active_level:
            if level_bed.metadata.cell_size < depth_h.level(depth_h.active_level).metadata.cell_size:
                initial = depth_h.prolong_to(level, _TRANSFER.depth_prolong_method)
            else:
                initial = depth_h.restrict_to(level, _TRANSFER.depth_restrict_method)
        else:
            initial = depth_h.level(level)
        initial = Raster.from_numpy(initial.to_numpy(), level_active, initial.metadata)

        options = dict(backend_options)
        options.update(
            point_inflows=_coarse_points(point_inflows, factor, bed.shape),
            initial_water_depth=initial,
            max_steps=max_steps,
        )
        result = solver(level_bed, level_state, level_rain, level_friction, **options)
        level_results.append(result)
        depth_h.set_active(result.water_depth)

    return GraphFloodMultigridResult(level_results[-1], tuple(level_results), factors)


_DR = np.asarray((-1, 0, 0, 1, -1, -1, 1, 1))
_DC = np.asarray((0, -1, 1, 0, -1, 1, -1, 1))


def _shift(values: np.ndarray, dr: int, dc: int, boundaries: str) -> np.ndarray:
    shifted = np.roll(values, shift=(-dr, -dc), axis=(0, 1))
    if boundaries[0] != "p":
        if dr < 0:
            shifted[: -dr, :] = 0
        elif dr > 0:
            shifted[-dr:, :] = 0
    if boundaries[1] != "p":
        if dc < 0:
            shifted[:, : -dc] = 0
        elif dc > 0:
            shifted[:, -dc:] = 0
    return shifted


def _dilate(mask: np.ndarray, cells: int, boundaries: str) -> np.ndarray:
    result = mask.copy()
    for _ in range(cells):
        grown = result.copy()
        for dr, dc in zip(_DR, _DC):
            grown |= _shift(result, int(dr), int(dc), boundaries)
        result = grown
    return result


def _outlet_connected(mask: np.ndarray, outlets: np.ndarray, boundaries: str) -> np.ndarray:
    rows, cols = mask.shape
    connected = np.zeros(mask.shape, dtype=bool)
    queue = deque(int(flat) for flat in np.flatnonzero(outlets & mask))
    for flat in queue:
        connected.flat[flat] = True
    while queue:
        flat = queue.popleft()
        row, col = divmod(flat, cols)
        for dr, dc in zip(_DR, _DC):
            nr, nc = row + int(dr), col + int(dc)
            if boundaries[0] == "p":
                nr %= rows
            if boundaries[1] == "p":
                nc %= cols
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if mask[nr, nc] and not connected[nr, nc]:
                connected[nr, nc] = True
                queue.append(nr * cols + nc)
    return connected


def _downstream_closure(
    mask: np.ndarray,
    surface: np.ndarray,
    original_active: np.ndarray,
    outlets: np.ndarray,
    boundaries: str,
) -> np.ndarray:
    result = mask.copy()
    while True:
        receivers = np.zeros(mask.shape, dtype=bool)
        for dr, dc in zip(_DR, _DC):
            neighbour_surface = _shift(surface, int(dr), int(dc), boundaries)
            neighbour_active = _shift(original_active, int(dr), int(dc), boundaries)
            donors = (
                result & ~outlets & neighbour_active &
                (neighbour_surface < surface)
            )
            receivers |= np.roll(
                donors, shift=(int(dr), int(dc)), axis=(0, 1)
            )
        expanded = result | (receivers & original_active)
        if np.array_equal(expanded, result):
            return result
        result = expanded


def _crossing_inflow(
    surface: np.ndarray,
    qin: np.ndarray,
    original_active: np.ndarray,
    induced_active: np.ndarray,
    cell_size: float,
    boundaries: str,
    exponent: float,
) -> np.ndarray:
    weights = []
    total = np.zeros(surface.shape, dtype=np.float64)
    for direction, (dr, dc) in enumerate(zip(_DR, _DC)):
        neighbour_surface = _shift(surface, int(dr), int(dc), boundaries)
        neighbour_active = _shift(original_active, int(dr), int(dc), boundaries)
        distance = cell_size if direction < 4 else cell_size * np.sqrt(2.0)
        downslope = original_active & neighbour_active & (neighbour_surface < surface)
        slope = np.where(downslope, (surface - neighbour_surface) / distance, 0.0)
        weight = np.where(downslope, np.power(slope, exponent) * distance, 0.0)
        weights.append(weight)
        total += weight

    crossing = np.zeros(surface.shape, dtype=np.float64)
    for weight, dr, dc in zip(weights, _DR, _DC):
        receiver_is_induced = _shift(induced_active, int(dr), int(dc), boundaries)
        donor_flux = np.divide(
            qin * weight,
            total,
            out=np.zeros_like(total),
            where=(total > 0) & ~induced_active & receiver_is_induced,
        )
        crossing += np.roll(donor_flux, shift=(int(dr), int(dc)), axis=(0, 1))
    return crossing


def _point_field(points: Any, shape: tuple[int, int]) -> np.ndarray:
    field = np.zeros(shape, dtype=np.float64)
    if points is None:
        return field
    for row, col, discharge in np.asarray(points, dtype=np.float64).reshape(-1, 3):
        field[int(row), int(col)] += discharge
    return field


def run_graphflood_induced(
    bed: Any,
    node_state: Any,
    rainfall_rate: Any,
    friction_coefficient: Any,
    qin_threshold: float,
    *,
    backend: str | Callable[..., Any] = "optimised_v2",
    halo_cells: int = 3,
    edge_depth_threshold: float = 1.0e-4,
    max_expansions: int = 2,
    point_inflows: Any = None,
    median_dhw_convergence_tolerance_metres: float | None = None,
    median_dhw_relative_convergence_tolerance: float | None = None,
    **backend_options: Any,
) -> GraphFloodInducedResult:
    """Run GraphFlood on a discharge-selected domain and floodplain halo.

    Parameters
    ----------
    bed, node_state, rainfall_rate, friction_coefficient : Raster
        Full-domain hydraulic inputs.
    qin_threshold : float
        Minimum full-domain volumetric ``Qin`` retained as channel domain.
    backend : str or callable, default="optimised_v2"
        Backend used by the prepass and reduced solves.
    halo_cells : int, default=3
        Initial dilation and wet-edge expansion distance.
    edge_depth_threshold : float, default=1e-4
        Depth that activates expansion at a retained-domain edge.
    max_expansions : int, default=2
        Maximum number of domain expansions after the first solve.
    point_inflows : array-like, optional
        Rows of ``(row, col, discharge_m3_s)``.
    **backend_options
        Options forwarded to the backend.

    Returns
    -------
    GraphFloodInducedResult
        Reduced result, prepass, active mask, and expansion history.

    Authors
    -------
    B. Gailleton
    """
    from . import Raster

    backend_options = _normalize_convergence_options(
        median_dhw_convergence_tolerance_metres,
        median_dhw_relative_convergence_tolerance,
        backend_options,
    )

    if qin_threshold < 0 or halo_cells < 0 or max_expansions < 0:
        raise ValueError("threshold, halo_cells, and max_expansions must be non-negative")
    solver, backend_name = _backend(backend)
    boundaries = backend_options.get("boundary_conditions", "oooo")
    exponent = float(backend_options.get("mfd_exponent", 1.0))
    prepass_options = dict(backend_options)
    prepass_options.update(max_steps=1, run_until_converged=False,
                           point_inflows=point_inflows)
    prepass = solver(
        bed, node_state, rainfall_rate, friction_coefficient, **prepass_options
    )

    states = node_state.to_numpy()
    original_active = bed.mask_numpy() & node_state.mask_numpy() & (states != 0)
    outlets = original_active & (states == 3)
    qin = prepass.qin.to_numpy().astype(np.float64)
    source_points = _point_field(point_inflows, tuple(bed.shape))
    # qin was routed on the surface before the prepass depth update.
    surface = (prepass.water_surface.to_numpy() - prepass.dhw.to_numpy()).astype(np.float64)
    channel = original_active & ((qin >= qin_threshold) | outlets | (source_points > 0))
    active = _dilate(channel, halo_cells, boundaries) & original_active
    active = _downstream_closure(
        active, surface, original_active, outlets, boundaries
    )
    active = _outlet_connected(active, outlets, boundaries)
    expansion_results = []
    initial_depth = None

    final_result = prepass
    for expansion in range(max_expansions + 1):
        crossing = _crossing_inflow(
            surface, qin, original_active, active,
            bed.metadata.cell_size, boundaries, exponent,
        )
        induced_states = np.where(outlets, 3.0, np.where(active, 1.0, 0.0))
        state_raster = Raster.from_numpy(
            induced_states.astype(np.float32), original_active, bed.metadata
        )
        rain_values = np.where(active, rainfall_rate.to_numpy(), 0.0)
        rain_raster = Raster.from_numpy(
            rain_values.astype(np.float32), bed.mask_numpy(), bed.metadata
        )
        point_field = np.where(active, source_points + crossing, 0.0)
        rows, cols = np.nonzero(point_field > 0)
        induced_points = np.column_stack((rows, cols, point_field[rows, cols]))

        options = dict(backend_options)
        options.update(point_inflows=induced_points,
                       initial_water_depth=initial_depth)
        final_result = solver(
            bed, state_raster, rain_raster, friction_coefficient, **options
        )
        expansion_results.append(final_result)
        if expansion == max_expansions:
            break

        outside_neighbour = np.zeros(active.shape, dtype=bool)
        for dr, dc in zip(_DR, _DC):
            outside_neighbour |= _shift(original_active & ~active,
                                        int(dr), int(dc), boundaries)
        wet_edge = active & outside_neighbour & (
            final_result.water_depth.to_numpy() > edge_depth_threshold
        )
        if not wet_edge.any():
            break
        expanded = _dilate(wet_edge, halo_cells, boundaries) | active
        expanded &= original_active
        expanded = _downstream_closure(
            expanded, surface, original_active, outlets, boundaries
        )
        expanded = _outlet_connected(expanded, outlets, boundaries)
        if np.array_equal(expanded, active):
            break
        previous_active = active
        active = expanded
        depth = np.where(
            previous_active & final_result.water_depth.mask_numpy(),
            final_result.water_depth.to_numpy(), 0.0,
        )
        initial_depth = Raster.from_numpy(
            depth.astype(np.float32), bed.mask_numpy(), bed.metadata
        )

    source_discharge = float(
        np.where(original_active, rainfall_rate.to_numpy(), 0.0).sum(
            dtype=np.float64
        ) * bed.metadata.pixel_area + source_points.sum(dtype=np.float64)
    )
    outlet_qin = float(final_result.qin.to_numpy()[outlets].sum(dtype=np.float64))
    relative_mass_error = abs(outlet_qin - source_discharge) / max(
        abs(source_discharge), 1.0e-30
    )
    return GraphFloodInducedResult(
        final_result,
        prepass,
        active,
        float(active.sum() / max(1, original_active.sum())),
        tuple(expansion_results),
        source_discharge,
        outlet_qin,
        relative_mass_error,
    )
