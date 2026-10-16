# Landscape Evolution Modelling in lsdtt3

This page explains the design choices available in lsdtt3's landscape evolution model (LEM): what incision and transport laws are implemented, what numerical schemes solve them, where the theory comes from, and what trade-offs guide you toward one choice or another. It is not a CLI reference — for parameter names and default values see [CLI reference: lsdtt-lem](../cli_reference/cli_reference_lsdtt_lem.md). Instead, this page is for understanding the landscape processes and numerics that each parameter controls.

It should be noted: lsdtt3's lem functionality has come about because it has been aimed at answering specific research questions and testing different components of the topographic analysis packages. It is not a full landscape evolution model with all the bells and whistles of, say, landlab. But if you just want a simulation using basic stream power + hillslopes on a rectangular grid lsdtt-lem will do the job. 

If you are modelling a landscape for the first time (I'm not sure you would be reading this page if that were the case, but whatever), start with [Understanding fluvial and hillslope processes](#fluvial-and-hillslope-laws-in-lsdtt3) and [Steady state versus time stepping](#steady-state-versus-time-stepping), then turn to the CLI reference to find parameter values that match your process choice.

---

## Overview

The lsdtt3 LEM computes landscape evolution by applying three processes to a DEM at each timestep (or iteratively for steady-state): **fluvial incision**, **hillslope sediment transport**, and **uniform tectonic uplift**. Each process is optional and can be switched on independently. The core idea is that erosion/deposition = stream power (on channels) + diffusive/threshold transport (on hillslopes) + uplift, and the landscape evolves by solving this differential equation forward in time or toward a steady state.

Two broad strategies are available:

- **Time stepping**: integrate the governing equations forward from an initial topography at small timesteps over long model duration.
- **Analytical steady state**: compute the fluvial part directly in one pass (no time loop), optionally solve the hillslope part iteratively, and optionally repeat ("snap") to let drainage divides migrate.

---

## Fluvial and hillslope laws in lsdtt3

### Stream power and fluvial incision

The **stream power law** describes erosion on channelized flow paths. It states that the rate of channel incision is

```
dz/dt = U − K · A^m · S^n
```

where `z` is elevation, `t` is time, `K` is an erodibility coefficient, `A` is the drainage area, `S` is the local slope, `m` and `n` are dimensionless exponents (typically 0.4–0.5 and 1, respectively, although n = 1 is more for computational convinience than for accurate physics, see the last two decades of landscape evolution papers), and `U` is uplift rate. Incision lowers the surface, so the erosion term enters
with a negative sign; the steady state `U = K A^m S^n` is unchanged by the
convention, but the transient is not. The area-slope form comes from
**Howard and Kerby (1983)**, who measured it directly on rapidly eroding badland channels. In lsdtt3 it is the basis of all fluvial incision.

#### The exponents m and n

- **`m`** (the drainage-area exponent, default 0.45) captures how incision scales with catchment size. A larger basin carries more water, so it erodes faster. `m ≈ 0.4–0.5` reflects observations that incision is sub-linear in area: doubling the catchment does not double the erosion rate. This exponent actually includes exponents about runoff efficiency as well as channel width scaling: see Whipple and Tucker (1999). 
- **`n`** (the slope exponent, default 1.0) determines how strongly the channel responds to slope. `n = 1` (linear slope dependence) is the simplest choice. This does have a physical basis (see Whipple and Tucker, 1999) but in reality it is often chosen due to its considerable numerical advantages over n not equal to 1. Empirical observations suggest n > 1 (e.g., Harel et al., 2016) but this is an ongoing debate within the geomorphology community, with many arguments for n = 1 and many against. `n ≠ 1` requires iterative (Newton–Raphson) solution instead of a closed-form formula, making it numerically (much) more expensive.

The combination matters for landscape form: the scaling of the characteristic
landscape length depends on the erosivity convention and on the characteristic
scales used to nondimensionalise the model. The two conventions therefore do
not share one universal landscape Péclet exponent.

#### Numerical schemes for stream power

The stream power law is stiff: steep hillslopes try to transport sediment very fast, and naïve time-stepping requires tiny timesteps. lsdtt3 solves this using **implicit (unconditionally stable) numerical schemes**:

- **n = 1 (default):** The implicit update is **linear** in the unknown and so has a closed form, `z_new = (z_old + K A^m (dt/dx) z_r) / (1 + K A^m (dt/dx))` — Braun & Willett (2013), eq. 2. This is the fastest route and is taken automatically when `n_exponent: 1.0`.

    One qualification: this applies to the **single-receiver** case. In the D-infinity path a cell with two receivers is solved by Newton–Raphson whatever `n` is, because the two-receiver equation is not linear even at `n = 1`.
- **n ≠ 1:** Newton–Raphson iteration over an implicit equation. More expensive but allows exponents like `n = 1.5` or `n = 2`. The implementation builds on the implicit stream-power solver strategy of **Braun and Willett (2013)** — the FastScape algorithm, an O(n) implicit sweep in topological order. **Anand, Hooshyar and Porporato (2020)** extend that ordering from single- to multiple-flow-direction networks, which is what makes a D-infinity fluvial step affordable at the same complexity; and **Anand et al. (2022)** give a linear stability analysis of when a fluvial landscape channelises and at what spacing, which is the analytical counterpart of the emergent drainage density discussed below.

See the CLI reference for the `n_exponent` and `slope_clamp` parameters that control this solver.

---

### Hillslope sediment transport

Hillslopes carry sediment differently from channels. The two regimes represented in lsdtt3 are:

#### Linear diffusive transport

The simplest model: sediment flux on a slope is proportional to gradient,

```
q = −D · ∇z
```

where `q` is the sediment flux (volume per unit width per time), `D` is a diffusivity coefficient (m² yr⁻¹), and `∇z` is the elevation gradient. This is a diffusion equation, and it is reversible — material spreads downhill and fills valleys, and at equilibrium the landscape is a parabola. It is named for **Culling (1960)**, who drew the analogy between hillslope erosion and Fourier heat conduction; **Culling (1963)** later derived the same form from random particle motion, which is the paper to cite when the mechanism rather than the equation is the point.

Advantages: simple, fast, always stable.
Disadvantage: real hillslopes often have a threshold slope beyond which transport accelerates sharply (e.g. soil creep at low slopes, landsliding at high slopes). Linear diffusion does not capture this.

#### Nonlinear (critical-slope) transport

Above a **critical slope** `S_c`, sediment transport becomes much faster — non-linearly so. The **Roering law** (Roering, Kirchner, and Dietrich 1999) captures this:

```
q = −D · ∇z / (1 − (|∇z| / S_c)^2)
```

At slopes much below `S_c`, this reduces to linear diffusion (the denominator → 1). As slope approaches `S_c`, the transport coefficient diverges, modeling an abrupt increase in sediment flux — a threshold effect. Above `S_c`, the model is often modified to ensure mass balance; a common approach is **Cidre-style** (Carretier et al. 2016, *Earth Surface Dynamics* 4:237–251), which applies a slope-restoration term to drive slopes that exceed `S_c` back toward `S_c`, so the landscape self-organizes toward a critical-slope topography.

This nonlinear law is physically motivated by hillslope processes: at low angles, soil creep dominates (nearly linear); at high angles, landsliding or soil raveling dominates (nonlinear). Real hillslopes rarely exceed `S_c ≈ 0.8` (45 degrees) in steady state.

---

## Numerical schemes and solvers for hillslope transport

lsdtt3 offers multiple solvers for hillslope transport. The primary choice is the **transient solver**, selected via `hillslope_model` (`cellular`, `cidre`, `roering_implicit`, or `adi`; legacy explicit names remain internal). For steady-state solves, Newton is used with SOR/Picard fallback.

> **Newton is the default steady-state solver.** Steady-state runs (`run_steady_state: true`)
> use Newton for both linear and nonlinear flux laws, and transient `roering_implicit` steps use
> the Newton backward-Euler stepper. No parameter is needed and `hillslope_sor_omega` no longer
> has to be tuned. SOR with Picard linearisation remains the automatic fallback — it is used when
> Newton is unavailable (a build without Eigen) or fails to converge, and can be forced with
> `hillslope_force_sor: true` for comparison against an older run.

**Transient solvers (selected via `hillslope_model`):**

| Solver | Equation | Time stepping | Stability | Best for |
|---|---|---|---|---|
| `linear` | D · ∇²z | explicit FD | CFL-limited | Thick sediment cover, long time horizons, when diffusion dominates |
| `roering` | Roering (1999) nonlinear law | explicit FD with auto-subcycling | CFL-limited near S_c | Transient dynamics, landscapes with strong slope feedback |
| `cidre` | Carretier et al. (2016) finite-volume | explicit FV + D-infinity routing | CFL-limited, but stable near Sc via slope-restoration | Coupled fluvial–hillslope on emergent network (D-infinity). Note that this does not have significant time savings over roering for test cases.|
| `roering_implicit` | Roering nonlinear law | implicit backward-Euler or TR-BDF2 (Newton, Picard fallback) | unconditionally stable (nonlinear) | Long, fine-grid runs; relief-targeting studies; see caveat below |
| `adi` | Laplacian | alternating-direction implicit | unconditionally stable (linear only) | Regular grids; when explicit schemes bottleneck |
| `cellular` | Roering below `r_match`, finite tangent continuation above | conservative face transfers with automatic subcycling | CFL controlled, finite above S_c | Fast nonlinear coupling, nonlocal transport and tracing |

### Linear diffusion solver (`linear`)

Solves `dz/dt = D · ∇²z` using an **explicit finite-difference stencil** (forward Euler). Stability is subject to the CFL condition: the timestep must not exceed `Δt ≤ h² / (4D)`, where `h` is the cell size. On fine grids (`h = 1` m) this can force very small timesteps unless `D` is large. However, it is simple, well-understood, and when the diffusivity is genuinely high (thick mantles of soil or saprolite), the CFL limit is not severe.

### Roering explicit solver (`roering`)

Solves the nonlinear Roering law using an **explicit finite-difference stencil** with **automatic subcycling** near critical slope. When a cell's slope approaches `S_c`, the time step is locally subdivided to remain stable. This is slower than linear diffusion but captures the nonlinear threshold behaviour without the cost of matrix solves. It is suitable for transient models where you want to see how hillslopes respond to a changing climate or base level, and where you have time to run.

### Cidre solver (`cidre`)

**Cidre** (from Carretier et al. 2016) is a **finite-volume forward sweep** along a D-infinity flow graph. It sweeps headwaters to outlet and routes sediment using receiver weights. It supports emergent drainage, but rebuilding/routing the graph can be expensive; `cellular` and transient `roering_implicit` also support an empty channel set.

It uses D-infinity routing (Tarboton 1997), so it is more expensive than a regular-grid stencil, but it scales well on distributed networks.

### Roering implicit solver (`roering_implicit`)

The **unconditionally stable backward-Euler solver** for the Roering law. It solves an implicit nonlinear equation with a sparse linear solver from Eigen3, and is the most robust for long runs on fine grids: being implicit, it carries no CFL restriction, so the timestep is set by accuracy rather than stability. (`adi` is also unconditionally stable, but only for **linear** diffusion; `roering_implicit` is the unconditionally stable option for the **nonlinear** law.)

**Newton first, Picard as fallback.** The nonlinear solve is Newton–Raphson by default: differentiating the nonlinear mobility gives an iteration count that does not grow with the grid. Picard is retained as the automatic fallback — it is the only option in an Eigen-free build, and a step Newton cannot take is worth retrying with a damped method rather than failing the run. `hillslope_force_sor` forces the Picard path for transient runs as well as steady ones.

Transient `roering_implicit` can run in the emergent-drainage regime: transient diffusion is well posed without channel Dirichlet nodes and distributed fluvial incision supplies the sink. Only a purely steady elliptic solve without any Dirichlet boundary is singular.

### Cellular soft-critical solver (`cellular`)

`cellular` computes one conservative flux per raster face. It follows Roering exactly below
`cellular_r_match` (default 0.85), then continues mobility with a finite tangent, allowing temporary
slopes above `S_c`. Near ridgetops it retains linear flux and the steady `U/D` curvature relation.
Optional nonlocal transport splits near-critical excess flux into a bounded geometric downhill
kernel; channels/open boundaries intercept it. Sparse Lagrangian superparticles sample the same
accepted fluxes without changing topography. Particle CSVs preserve provenance, age, distance,
random seed and substep index for deterministic restart.

### ADI solver (`adi`)

**Alternating-Direction Implicit** (ADI) method for solving the Laplacian, splitting the 2D problem into two 1D line solves and repeating. It is unconditionally stable and very fast on regular grids, but it applies only to linear diffusion (no nonlinearity). It does not require Eigen and is a useful fallback when explicit schemes bottleneck but you want to avoid a full sparse-matrix solve.

### The steady-state hillslope solver: red–black SOR

The five solvers above all march forward in time. When `run_steady_state: true` the hillslope
field is instead obtained directly, by solving the **steady** equation in which uplift is
exactly balanced by divergence of the sediment flux:

    ∇·q_s = U          with    q_s = -D ∇z / (1 - |∇z|²/S_c²)

with channels held as Dirichlet boundaries. Discretised on the regular grid this is a
system of `N = n_rows × n_cols` coupled equations, and it is solved iteratively by
**red–black successive over-relaxation (SOR)**.

#### What SOR is, and why over-relaxation

Take the linear case (`S_c → ∞`, so the flux law reduces to `D ∇²z = -U`). The five-point
stencil makes each cell's elevation the average of its four neighbours plus a source term:

    z_jac = (z_N + z_S + z_E + z_W + U h²/D) / 4

Sweeping that assignment over the grid is **Gauss–Seidel** relaxation. It converges, but
slowly: each sweep propagates information only one cell, so long-wavelength error decays at
a rate set by the domain size. **Over-relaxation** deliberately overshoots the Gauss–Seidel
update by a factor `ω`:

    z_new = (1 - ω) z_old + ω z_jac

For `ω = 1` this is Gauss–Seidel; for `1 < ω < 2` it extrapolates past the new estimate,
which damps the slowest error modes far more aggressively. The optimal value follows from
the Jacobi spectral radius `ρ_J`:

    ρ_J = ½ [cos(π/(N+1)) + cos(π/(M+1))]       ω_opt = 2 / (1 + √(1 - ρ_J²))

`ω > 2` diverges unconditionally. lsdtt3 computes `ω_opt` from the grid dimensions unless
you override it with `hillslope_sor_omega`.

**"Red–black"** refers to the ordering. Colour the grid like a chessboard: a cell's four
von Neumann neighbours are always the opposite colour. So every red cell can be updated
using only black values, and vice versa. Each iteration is therefore two colour passes, and
within a pass every update is independent — which is what makes the sweep safe to
parallelise with OpenMP, unlike a naive lexicographic Gauss–Seidel sweep.

#### The expected cost: O(n) iterations

This is the property to check a run against. For Gauss–Seidel the error per sweep decays by
`ρ_GS ≈ 1 - O(1/n²)` on an `n × n` grid, giving **O(n²) iterations** to a fixed tolerance.
Optimal SOR improves the decay rate to `ρ_SOR ≈ 1 - O(1/n)`, giving **O(n) iterations** —
the square root of the Gauss–Seidel count, and the whole reason for using it.

Concretely, for `n = 250`: `ρ_SOR ≈ 0.975`, so reaching a tolerance of `10⁻⁴` needs roughly
`log(10⁻⁴)/log(0.975) ≈ 364` iterations. A 250 × 250 steady-state solve that takes a few
hundred iterations and a couple of seconds is behaving correctly. **One that runs for tens
of thousands of iterations is reporting a defect, not doing hard work** — either the
relaxation factor is wrong, the problem is not the one SOR is being asked to solve, or the
convergence test itself is broken (below).

#### Picard iteration for the nonlinear law

With a finite `S_c` the problem is no longer linear: the effective diffusivity at each cell
face depends on the slope there, which depends on the very elevations being solved for. The
face coefficient is

    α_face = 1 / (1 - |∇z|²_face / S_c²)

which rises steeply as the slope approaches `S_c`. lsdtt3 handles this by **Picard
iteration** (also called successive substitution, or lagged-coefficient iteration): at each
sweep the coefficients `α` are evaluated from the *current* elevation field, then held fixed
while the cell is updated as if the problem were linear. The updated field then supplies new
coefficients on the next sweep. Repeat until the field stops changing.

Picard is attractive because it reuses the linear machinery unchanged and needs no Jacobian
(unlike Newton's method, which converges quadratically but requires differentiating `α` with
respect to every neighbouring elevation). Its weakness is that it has **no damping on the
coefficients**. Under-relaxing `ω` damps the change in `z`, but `α` is still recomputed at
full strength from whatever `z` currently is. On steep faces, where `dα/dS` is large, this
permits a feedback loop: the elevation drops, the face steepens, `α` rises, the update pulls
harder, and the cell overshoots — then the same thing happens in reverse. The result is a
**limit cycle** rather than convergence, confined to the small fraction of cells sitting at
high `S/S_c`. Nonlinear runs that report a residual which stops falling and simply hovers,
while the *field* is visibly steady, are showing this signature.

#### Why this solver needs double precision

The SOR iteration is carried in **double precision** even though the rasters either side of
it are single precision. This is not incidental tidiness; the solver does not work without
it, and the reason generalises to any iterative solver judged on an increment.

The convergence test is `max |z_new - z_old| < tolerance`: a quantity that by construction
shrinks toward zero, evaluated on a field that does not. Elevations are O(100 m), and the
spacing between representable `float` values near 100 m is about `10⁻⁵` m. Worse, the update
forms a sum of four neighbours — O(400) in magnitude, where the spacing is ~`6 × 10⁻⁵` m —
before differencing it back down against `z_old`. The subtraction is catastrophically
cancelling: nearly every significant digit of the operands is common to both.

So in single precision the measured increment cannot fall below roughly `2 × 10⁻⁴` m no
matter how well converged the solution is. It bottoms out there and then **wanders**,
driven by rounding rather than by the physics. Since the default tolerance is `10⁻⁴` m, the
test sat permanently below the arithmetic noise floor: `converged` could never become true,
and every solve ran to `hillslope_ss_max_iterations` regardless of its actual state. On a
250 × 250 grid this turned a 261-iteration, ~2-second solve into a 250,000-iteration,
13-minute one that still reported failure. The residual is not merely imprecise in this
regime — it is *non-monotonic*, so it cannot be used to compare two runs either.

Double precision moves that floor down by about nine orders of magnitude, which puts it far
below any tolerance a user would set and restores the increment's meaning in metres. The
memory cost is one extra array of `N` doubles, which is negligible against the iteration
count it saves. The elevation raster is converted back to `float` on the way out; only the
iteration is carried at the higher precision.

The general lesson: **choose the precision of an iterative solve from the smallest quantity
you intend to measure, not from the size of the field.** A tolerance is only meaningful if
it sits above the arithmetic noise floor of the expression that computes it, and a
convergence criterion built on a difference of large, nearly-equal numbers has a much higher
floor than the storage format alone would suggest.

### The shared discrete operator

All hillslope solvers in lsdtt3 — SOR steady-state, explicit transport, ADI, implicit transient, and Newton — operate on the same **conservative cell-centred finite-volume discretisation** on a regular grid. For an interior hillslope cell P with four cardinal neighbours, the governing equation is

```
dz_P/dt = U + (1/h²) Σ_f D_eff,f (z_nb,f - z_P)
```

where U is uplift rate, h is the cell size, and the sum runs over the four faces (N, S, E, W). The effective diffusivity at each face is

```
D_eff,f = D / (1 - (S_n,f² + S_t,f²)/Sc²)
```

where S_n is the slope component normal to the face, S_t is the slope component tangential to the face, and S_c is the critical slope. In the linear case (S_c → ∞), this reduces to the standard Laplacian with scalar diffusivity D.

**Critical point: face symmetry.** The face diffusivity is a property of the *face itself*, not of the cell being updated. This means the tangential slope must be **averaged across the face** from the two cells sharing it. If instead you evaluate it at only the cell being updated — say, using the E–W cross-slope at cell P for both the N face of P and the S face of its northern neighbour — the operator becomes non-symmetric: flux leaving P through a face does not equal flux entering the neighbour through it. Non-symmetric operators admit complex eigenvalue pairs, leading to oscillating rather than decaying modes, and the scheme becomes non-conservative. Measurements (see [testing_and_benchmarking.md §10e](../developer_notes/testing_and_benchmarking.md)) show this defect caused a **575% flux imbalance** in the old asymmetric code. The correct face-symmetric form is now shared by all solvers.

### Boundary conditions for hillslope solvers

Boundary conditions are specified as a four-character string (`boundary_conditions` parameter), one character per cardinal direction (N, E, S, W in that order). A short string pads with `'o'` on the right. Each character takes one of three values:

- **`'p'` — Periodic wrap.** The cell wraps to the opposite side, so a cell at the eastern edge has its eastern neighbour at the western edge. Used for cyclic domains to eliminate edge effects.
- **`'o'` — Dirichlet (fixed elevation at base level).** The edge cell is held at its input elevation and is not updated; it acts as a fixed-elevation boundary condition. Channel cells are always treated as `'o'` boundaries regardless of the boundary string.
- **`'n'` — Neumann zero-flux (reflecting).** The edge cell is reflected: the gradient normal to the boundary is set to zero by treating the edge as a mirror, so sediment flux across the boundary is zero.

The default is `"oooo"` (all edges fixed at base level), which is appropriate for a confined study area. A periodic domain uses `"pppp"` or (e.g.) `"poop"` for periodic east–west with fixed north–south edges.

### Newton's method for nonlinear steady-state

For the nonlinear Roering law with a finite critical slope, the nonlinear SOR solver (Picard iteration) has a fundamental limitation: **it has no damping on the coefficients**. Under-relaxing the update (reducing ω) damps the change in elevation `z`, but the diffusivity coefficient α is still recomputed at full strength from the current `z`. On steep faces where `dα/dS = α²` is large, this allows a feedback loop: elevation drops → face steepens → α rises → update pulls harder → overshoot. The result is a **limit cycle** confined to cells near the critical slope, and Picard cannot escape it at any fixed relaxation factor.

A nonlinear steady-state solve using Picard at optimal SOR relaxation (`ω ≈ 1.95`) or higher diverges; at `ω = 0.3` it converges in ~5,100 sweeps but is a 500× penalty against the linear optimal-SOR performance. Measurements (see [testing_and_benchmarking.md §10f](../developer_notes/testing_and_benchmarking.md)) on a 41×41 cross-slope landscape show SOR at `ω = 1.5` failing (residual plateaus at 26.6 m), whereas **Newton converges in 42 iterations to `||F||_∞ = 4e-17` m/yr** (an erosion-rate imbalance at machine precision).

**The core insight:** the steady problem and the backward-Euler transient step are the same root-find problem with one term:

```
F_P(z) = (z_P - z_P^n)/Δt - U - (D/h²) Σ_f α_f(z)(z_nb,f - z_P) = 0
```

Setting `1/Δt = 0` recovers steady-state; non-zero `1/Δt` is the transient form. A single **Newton core**, parameterised by the time-step reciprocal, handles both. This means the steady solver and the transient integrator use identical numerics and cannot drift.

**The Jacobian** of this system is 9-point (the cross-slope terms in `α_f` reach the diagonal neighbours) and non-symmetric (the `∂α` terms multiply `(z_nb - z_P)`, which is antisymmetric). This rules out the symmetric SPD solver (`SimplicialLDLT`); lsdtt3 uses `SparseLU` for moderate grids or `BiCGSTAB` with `IncompleteLUT` preconditioning for large systems.

**Convergence criterion:** Newton is judged on the residual **`||F||_∞`** in units of m/yr — a physically interpretable erosion-rate imbalance, independent of elevation magnitude. This avoids the catastrophic-cancellation trap that plagued the SOR solver: a criterion built on `max|Δz|` has a very high noise floor when the differences are computed from sums of O(100)-magnitude neighbours. Report both the residual norm and the max elevation change; a solve that looks good on one but not the other is hiding information.

**Globalisation and feasibility:** a pure Newton step can overshoot the critical slope. Two safeguards are applied:

1. **Line search on the residual norm:** if `||F(z + s)|| > ||F(z)||`, backtrack the step `s ← s/2` to satisfy descent. This ensures monotonic residual reduction.
2. **Feasibility clip:** before testing descent, check whether the trial step would put any face slope above `r_max·Sc` (where `r_max ≈ 0.99`). If yes, shrink the step until it does not. This prevents Newton from stepping over the critical slope barrier. A subtlety: if the initial field is already super-critical, a strict "never leave the feasible region" rule deadlocks (backtracking drives the trial toward that same infeasible field). The filter is conditional: from a feasible state it refuses to leave; from an infeasible one it accepts any restoration step that reduces the slope ratio, until the field crosses back below Sc.

**When the cap activates:** a run with `r_max` binding is a run where behaviour near Sc is determined by the limiter, not the flux law. This region is physically important (it is where critical-slope hillslopes regulate their form), so any run should flag and report where this occurs. On measured landscapes it is rare; the OCN test landscape reaches `max(S/Sc) = 0.882`, leaving the `r_max = 0.99` cap inactive everywhere.

**Newton is faster than SOR for the LINEAR problem too, which is not obvious.** A linear steady
solve has no critical-slope difficulty at all, so the natural expectation is that SOR wins. It does
not: `F` is already linear in `z`, so the Jacobian *is* the operator and Newton reduces to a single
sparse direct solve. The measured grid-by-grid comparison is in
[testing_and_benchmarking.md §10h](../developer_notes/testing_and_benchmarking.md),
which is the durable record for it; the numbers are not duplicated here. The
direction of the result is that Newton's advantage **grows with grid size**,
from near-parity on the smallest grid to more than an order of magnitude on the
largest tested.

Two things are going on. First, SOR's sweep count grows as roughly O(n²) here, not the O(n) that
over-relaxation is supposed to deliver, because the automatic `ω` is capped at 1.95 while the
optimum rises past it (1.962 at n = 161, 1.985 at n = 401). Above about n = 160 the cap binds and
SOR degrades toward Gauss–Seidel. Second, and more subtly, the two solvers are not delivering the
same accuracy: Newton's direct solve is exact to round-off, and the gap between them grows to
9.8e-4 m at n = 401 despite SOR being asked for 1e-6. That gap *is* SOR's error — an iterative
step-size criterion understates the true solution error by the usual `1/(1 − ρ)` factor, which is
about 1000× here.

**Where a direct solve stops paying.** Newton's cost is dominated by sparse-LU fill-in, and the
time exponent in `N` sits near 1.5 up to a few hundred thousand unknowns before climbing (the
measured points are in the benchmarking record). Above `direct_solve_max_unknowns` (default
250,000) the solver therefore switches automatically to BiCGSTAB with an incomplete-LU
preconditioner, which trades the factorisation for iterations and keeps memory bounded.

**Measurements (41×41 cross-slope fixture):** outer-iteration counts of 10 / 42 / 6 for n = 21 / 41 / 81 — the point being that the count does **not grow with the grid**, unlike SOR's O(n). It is not monotonic in n either: the count depends on how far each initial guess sits above `S_c`, because the restoration phase pays for that (the n = 41 case spends 89 line-search backtracks getting back under `S_c`). Linear mode solves in exactly one iteration and matches SOR to 3.8e-6 m; the transient form relaxes to the steady state in 15 sub-steps (49 Newton iterations) to within 1.9e-6 m.

### Stability: explicit versus implicit schemes

**Explicit (forward Euler) schemes are CFL-limited.** For the linear Laplacian, the stability condition is `Δt ≤ h²/(4D)`. For the nonlinear Roering law, the limit depends on the **tangent modulus** `dq/dS = D(1+r²)/(1-r²)²`, not the secant effective diffusivity `D_eff = D/(1-r²)`. This distinction matters: on a slope at `r = 0.93` (close to critical), the tangent modulus is ~13× larger than the secant, so the true CFL limit is ~13× more stringent than using `D_eff` in the formula would suggest. A solver that computes the CFL limit from `D_eff` rather than the tangent modulus reports a maximum stable timestep far larger than the scheme can actually tolerate.

The explicit Roering solver (`roering`) handles this by **auto-subcycling**: cells near `Sc` are sub-stepped during the global timestep to remain stable. This is reliable but computationally expensive on fine grids where many cells are steep.

**Implicit (backward-Euler) schemes are unconditionally stable.** There is **no CFL condition**. The timestep is chosen for *accuracy* and *nonlinear solvability*, not for stability. In a transient Newton solver with adaptive sub-stepping, a rejected sub-step means the nonlinear system was too hard to solve at that size or an accuracy gate fired — it never means the scheme went unstable. This is why implicit schemes are strongly preferred for long, fine-grid runs: they trade dense matrix solves for the freedom to take large timesteps without numerical instability. The cost-benefit depends on the problem; on short runs or coarse grids explicit schemes may be faster, but on geological timescales or sub-metre resolution the unconditional stability of implicit schemes is essential.

**Practical guidance:** choose `cellular` for fast nonlinear transient runs, emergent drainage, nonlocal transport or tracing. Choose `roering_implicit` when exact Roering behavior near its subcritical barrier is required and sparse-solve cost is acceptable. Newton remains recommended for nonlinear steady-state work.

---

## Comparison table: all hillslope solvers

Here is a comprehensive view of all available options, with the trade-offs that guide solver choice:

| Solver | Equation | When to use | Stability | Cost model | Notes |
|---|---|---|---|---|---|
| `linear` (explicit) | D·∇²z | Transient, thick diffusivities, if CFL is not limiting | CFL: Δt ≤ h²/(4D) | ~1 ms per cell per step | Simplest, well-tested; use if explicit schemes are fast enough |
| `roering` (explicit) | Roering with auto-subcycling | Transient, nonlinear dynamics, modest grids | CFL: Δt ≤ h²/(4D_tangent) near Sc | ~10 ms per cell per step; auto-subcycling adds cost near Sc | Captures nonlinear threshold; no matrix solves; expensive on fine grids |
| `cidre` (explicit FV) | Roering on D-infinity flow net | Emergent drainage density, coupled fluvial–hillslope | CFL-limited | ~100 ms per sweep (D-infinity routing cost) | Most realistic for channelization studies; used with an empty channel set |
| `adi` (implicit) | Laplacian only | Long linear runs on regular grids | Unconditionally stable | ~10 ms per line solve (1D subproblems) | Very fast; no sparse matrix fill; linear only |
| `roering_implicit` (backward-Euler, Newton with Picard fallback) | Roering nonlinear law | Long runs with imposed channels; relief targeting | Unconditionally stable | ~20–40 s per steady solve (62.5k unknowns); transient uses similar per step | Newton needs no relaxation factor; the Picard fallback is limited to ω ≤ 0.5 on nonlinear and limit-cycles at high ω |
| **Newton** (backward-Euler + Newton) *(default for steady state)* | Roering nonlinear law | Nonlinear steady-state and transient, especially near critical slopes | Unconditionally stable | ~10–30 s per steady solve (62.5k unknowns); O(1) iterations independent of grid | Converges where Picard cannot; 9-point Jacobian costs more per iteration but far fewer iterations; use SparseLU for grids <10⁶, BiCGSTAB/ILUT for larger |
| `cellular` | Soft-critical conservative flux, optional nonlocal split | Fast transient coupling and tracing | Explicit tangent-based subcycling | 19–136x faster than Newton in tested one-step grids | Exact Roering below `r_match`; finite above; deterministic topography |
| SOR steady-state (Picard) | Roering nonlinear law | Linear steady-state, or nonlinear at low ω (legacy) | N/A (steady-state only) | ~5–10 s per solve; O(n) iteration count scales with grid | Same operator as above; serves as fallback/predictor for Newton; documented here for completeness |

**Choosing a steady-state approach:**

- **Linear hillslopes:** use SOR (`run_steady_state: true`, `enable_nonlinear_hillslopes: false`). It is fast and exact.
- **Nonlinear hillslopes, imposed channels:** `run_steady_state: true`, `enable_nonlinear_hillslopes: true`. Newton is used automatically and needs no relaxation factor.
- **Legacy code or special cases:** the Picard solver within `roering_implicit` is still available as a fallback, but it should not be a first choice for nonlinear work.

---

## Steady state versus time stepping

lsdtt3 offers two fundamentally different routes to evolve a landscape: **transient (time-stepping)** and **analytical steady-state with optional iteration**.

### Time-stepping mode (`run_steady_state: false`)

The straightforward approach: apply fluvial incision, hillslope transport, and uplift as a sequence at each timestep, advance by `timestep` (yr), and repeat for `number_of_iterations` steps or until `stop_at_steady_state` is met.

This is the mode that directly solves the governing differential equations and is most appropriate when you want to **see transient dynamics** — e.g., a landscape's response to a change in climate, tectonics, or base level. It is also the mode that can run the coupled emergent-drainage case with Cidre (or other explicit schemes).

A landscape reaches **statistical steady state** (mean erosion rate = uplift rate) after a time that depends on the basin size and process rates. For a small synthetic model, this might be a few hundred timesteps; for a large real landscape, millions of years.

### Steady-state mode (`run_steady_state: true`)

Rather than integrating forward in time, the code uses an **analytical chi inversion** to compute the steady-state fluvial elevation profile in one pass:

1. **Integrate chi** along every D8 flow path, starting from the outlet and going upslope, where chi is `∫ A^(−m/n) ds` (an integral along the path). The steady-state elevation is obtained by inverting the stream power law.

2. **Optionally iterate** (the "snapping cycle"): take the computed elevations, recompute the flow graph (so divides can migrate), and snap again. This is necessary because a single pass leaves sharp steps across divides — each divide side is built independently from its own channel, with nothing forcing the two sides to meet smoothly.

3. **Optionally solve hillslope steady state** after fluvial snapping: on the final fluvial surface, set channel pixels as fixed Dirichlet boundary conditions and solve the hillslope equation iteratively to reach steady state (solving `D ∇²z = U` for linear, or the Roering equation for nonlinear).

#### Why snapping: divide equilibration

A **single-pass snap** leaves an elevated step along drainage divides. This is not an error — it is a fundamental consequence of building each side of a divide independently. The solution is to **iterate**: after snapping, recompute flow routing from the new elevations (so the divides can migrate if they want to), and snap again. Over several cycles, the divides settle into an equilibrium position where the two sides agree.

You can control the iteration using `steady_state_n_cycles` (how many to run) and optionally `steady_state_check_convergence` (stop early if the landscape stops changing). The equilibration cycle (enabled by `use_landscape_equilibration_cycle: true`) adds detailed diagnostics — per-cycle divide metrics, per-pixel divide steps, and per-cycle routing changes — so you can inspect whether the landscape actually settled.

**Caution:** A landscape need not converge to a fixed point. It can settle into a repeating **limit cycle** (e.g., oscillating drainage configuration) or wander indefinitely (especially with critical-slope hillslopes). The diagnostics warn you if a cycle repeats; if the reader's eye can't see a fixed point after many cycles, you may not have one.

#### Relief-mode calibration

One special steady-state feature: `steady_state_mode: relief`. Instead of specifying `erodibility` (K), you specify `target_relief` and the code back-calculates the K that produces that relief. This is convenient for target-based studies — "I want a landscape with relief between X and Y metres" — without having to guess K.

**Caveat:** The maximum relief is computed as `max_chi_on_domain` (over all cells), so on a **fluvial-only landscape** (no hillslopes), the channel network extends all the way to the divides. Channel tips are inherently steeper than real hillslopes, so the maximum relief is higher than it would be with hillslopes engaged. **Introducing hillslope transport lowers the relief below the target**, because the hillslopes relax the steep channel tips. If you use relief-mode with intended hillslope coupling, iterate: solve fluvial steady state with relief-mode, then solve hillslope steady state, measure the relief, and re-run fluvial steady state with a higher target.

---

## Spatially varying erodibility

Rock strength varies in space, so `lsdtt-lem` (`erodibility_fname`) and the
Python `LandscapeModel` (`erodibility_field`) accept a per-pixel erodibility
K(x). The field holds **absolute** K in the run's area convention; the scalar
`erodibility` stays as the reference K that a forcing history scales relative
to (K(x, t) = K_field(x) K_forcing(t) / K_ref).

### In the implicit fluvial step

The implicit solvers (D8 and D-infinity, closed form at n = 1 and Newton–Raphson
otherwise) use K in exactly one place, the coefficient K·a^m of each pixel's
update. With a field, pixel i uses K_i there; every update path and both
overexcavation clamps read that one coefficient. A uniform field equal to the
scalar therefore performs the same arithmetic on the same values, so its result
is bitwise identical to the scalar run, which the tests assert for D8 and
D-infinity at n = 1 and 2. The field is read afresh on every step, so a caller
may rewrite it between steps (for example from the elevation of each pixel, the
mechanism a lithology model needs); it is evaluated at the start of the step.

### In the analytical steady state: K-weighted chi

At steady state U = K_i a^m S^n on every pixel. Writing K_ref for the scalar,

```
z_i    = z_b + (U / K_ref)^(1/n) · χ_K,i
χ_K,i  = χ_K,recv + dx_eff · w_i · a_i^(−m/n),    w_i = (K_ref / K_i)^(1/n)
```

so a harder pixel (small K_i, large w_i) adds more chi and stands higher. In
D-infinity cells with two receivers the same weight multiplies a^(−m/n) both in
the discriminant and in the cardinal fallback. The solution does not depend on
K_ref (it cancels), and with w_i = 1 the arithmetic is the scalar one exactly.
In relief mode the field is scaled by one factor, `erodibility_field_scale`
= K_used / K_ref, so the relief reaches the target while the pattern of K is
kept.

### What the tests establish

* **Vertical contacts.** Along a channel crossing a contact from K1 into K2, the
  steady slope satisfies S2/S1 = (K1/K2)^(1/n) (Perne et al. 2017, their Eq. 4,
  for contacts that do not move with the rock). On an 80 × 3 strip with m = 0,
  U = 1e-3 m/yr and K2 = K1/4, time-stepping to steady state gives
  S1 = 0.1000 and S2 = 0.39999 against the expected 0.1 and 0.4 at n = 1, and
  0.316225 / 0.632452 against 0.316228 / 0.632456 at n = 2.
* **Every pixel.** On red-noise surfaces with a two-region field, run until no
  pixel moves more than 1e-3 U dt per step, K_i a^m S^n = U holds on every pixel
  to within 1e-3, for D8 and D-infinity, n = 1 and 2. The identity must be
  measured on the surface the implicit step solved (before uplift and before the
  boundary pixels are restored): boundary pixels that are not flow outlets erode
  within the sweep and are reset only afterwards, so their donors are solved
  against the eroded value.
* **The analytical steady state.** On a 200-node chain with a fourfold K
  contrast the K-weighted chi reproduces a double-precision reference to within
  200 ulp of the maximum elevation, and U = K_i a^m S^n holds link by link to
  about 5e-6.

**The snapped surface is not a fixed point of the transient model.** The
analytical snap integrates chi on the routing of the surface it is given, and
the snapped surface routes differently: on a 30 × 30 test surface about 330 of
840 D8 receivers change, with a uniform K as much as with a field. Time-stepping
the snapped surface then settles into a limit cycle (the largest per-step change
stays at 0.015–0.04 m) rather than a strict steady state. A field snap is as
close to steady as a scalar one (the first step after snapping moves the
surface by the same 99.9th-percentile amount, 1.368 m against 1.365 m); a region
of low K then takes proportionally longer to settle.

---

## Erosivity conventions: Perron versus Bonetti

A seemingly innocuous parameter, `divide_area_by_dx`, actually encodes a choice between two different conventions for what "drainage area" means in the stream power law. **This choice affects both the interpretation of `incision_threshold_specific_area` and the scaling of model results.** Getting it wrong is an error in an exponent, not a fudge factor, and cannot be calibrated away.

### The two conventions

**Convention 1 (Bonetti et al. 2020): Specific area** (`divide_area_by_dx: true`, default)

Compute `a = A / dx` (drainage area divided by cell size), so `a` has units of **metres**. The stream power law becomes

```
dz/dt = U − K · (a)^m · S^n  =  U − K · (A/dx)^m · S^n
```

This is the **default in lsdtt3** and was advocated by **Bonetti, Hooshyar, Camporeale and Porporato (2020, PNAS 117:1375–1382)** and **Litwin et al. (2022, JGR-ES)**.

Note that the units of `K` are `m^(1 − m) yr^-1`, which differs from
Convention 2. Because a channel's drainage area A is physical, a specific-area K
is tied to the cell size: `K_specific = K_map · dx^m`. A K calibrated on a 10 m
grid is not the same erodibility on a 30 m grid; `lsdtt-lem` converts a K raster
whose calibration cell size or convention differs from the run's, and prints the
factor.

**Convention 2 (Perron et al. 2008, 2009): Total area** (`divide_area_by_dx: false`)

Use the raw upslope area `A` in m², without dividing by dx. The law is

```
dz/dt = U − K · A^m · S^n
```

This was used in **Perron et al. (2008, JGR-ES 113:F04016)** and **Perron et al. (2009, Nature 460:502–505)**.

Note that the units of `K` are `m^(1 − 2m) yr^-1`, which differs from
Convention 1.

### Two critical consequences

**1. Units of the incision threshold**

- **`divide_area_by_dx: true`** → `incision_threshold_specific_area` is **A/dx in metres**. A value of 1,000,000 m means "cells with at least 1 million metres of upstream area per cell width". On a 10 m DEM this is 100,000 pixels; on a 1 m DEM this is 1,000,000 pixels.
- **`divide_area_by_dx: false`** → `incision_threshold_specific_area` is a **map area in m²**. A value of 1,000,000 m² means 1 km², regardless of cell size.

**2. Scaling of model results**

The characteristic length and height scales of the model depend on the
erosivity convention. For the linear hillslope-diffusion model, the relevant
scalings are:

- **Specific area** (Bonetti convention): use the `L_g` formulation of Litwin
  et al. (2025), with `a = A/dx` and `K` interpreted using the specific-area
  convention.
- **Total area** (Perron convention): use the `L_h` formulation of Theodoratos
  et al. (2018), with raw area `A` and `K` interpreted using the total-area
  convention.

Thus, the scaling of model results is different under the two conventions. A
landscape designed with the same numerical values of `m`, `n`, `K`, and `D`
under the two conventions does not in general have the same characteristic
length or height scale. The convention, the coefficient units, and the chosen
`dx` must be stated when comparing model results.

These scaling relationships apply only to models with **linear hillslope
diffusion**. Nonlinear or threshold hillslope laws introduce additional scales
and cannot be interpreted using these linear-diffusion characteristic-scale
relations alone. Although such models can be run with `lsdtt-lem`,
`lsdtt-advection-diffusion` should preferably be used for clean numerical
experiments with advection–diffusion models intended to follow these scaling
relationships. It provides a deliberately minimal fixed-timestep model with
linear ADI diffusion, controlled fluvial incision, and uniform uplift.

### Choosing a convention

If you are setting up a new model, **use the default (`divide_area_by_dx: true`, specific area)**, which matches recent literature and is the lsdtt3 convention. If you are reproducing a published study or comparing against code that uses raw area, flip `divide_area_by_dx: false` and adjust your `incision_threshold_specific_area` and any calibrated constants accordingly.

The CLI reference documents the resolved threshold value in the run header so you can double-check what was actually compared.

---

## Imposed versus emergent drainage density

One of the most important design choices in landscape evolution modelling is whether the **channel network is fixed or allowed to emerge**. lsdtt3 supports both.

**Say which one you mean.** `hillslope_channel_coupling` takes `auto`, `coupled` or
`channels_as_boundaries`. Historically the choice was an unreported side effect of
`incision_threshold_specific_area`, and a run could not be interpreted afterwards without
re-reading its parameter file. `auto` is the default and reproduces that inference
exactly — a positive threshold means channels are boundaries — but it now *reports*
what it resolved to, and the explicit values override the threshold rather than
merely describing it.

These are two different experiments, not two settings, so the combinations that
contradict themselves are refused rather than silently resolved: asking for
`coupled` with `enable_d8_fluvial_erosion: true` or with `run_steady_state`;
asking for `channels_as_boundaries` with no channel set at all; and asking for
`roering_implicit` with `run_steady_state` when `boundary_conditions` contains
no `'o'` edge — that last case is a steady elliptic problem with no Dirichlet
data anywhere, so the solution is not determined and the field drifts without
bound.

`coupled` together with `roering_implicit` **is** allowed. It used to be fatal;
that restriction was removed on 2026-08-27 (see below).

### Imposed drainage density (`incision_threshold_specific_area > 0`)

When you set `incision_threshold_specific_area > 0`, the code identifies channels as any cell meeting the resolved A/dx or map-area threshold. Channels become **Dirichlet boundary conditions** for the hillslope equation — they have fixed elevation during the hillslope solve and act as sediment sinks. The drainage density is an **input**, not an output: you choose the network spacing by picking the threshold.

Mind the units when picking a value — they depend on `divide_area_by_dx` (see
above). On a 10 m grid, a 10 km² network is `10 000 000` in map area, or
`1 000 000` as specific area; the two numbers differ by a factor of `dx`, and the
run header prints all three forms so a log can be read without the parameter file.

This is the approach of the **MUDDPILE model** (from which lsdtt3 was partly descended) and is appropriate when you have a **field-mapped channel network** or when you are studying hillslope form **conditional on a known network**.

The channel set is computed once per fluvial solve (D8 routing) and remains fixed during that timestep's hillslope transport. The network is fast and deterministic.

### Emergent drainage density (`incision_threshold_specific_area = 0`)

When `incision_threshold_specific_area: 0`, there are no pre-defined channels: **both fluvial and hillslope processes act everywhere**, and the network must **emerge**. Sediment is partitioned between fluvial and diffusive transport based on local slope and area, and the landscape self-organizes into channels where flow concentrates.

This is the approach of the **DonkModel** and is physically more realistic — channel heads are not imposed but arise where hydrology and sediment dynamics favor their development. However, it is **more sensitive to numerics and requires care**:

- **D-infinity routing is required** for hillslope coupling to work correctly. With D8, a diverging hillslope cell (upslope area increases downslope) may manufacture spurious channel heads because area concentrates artificially along one receiver.
- **The coupled regime is slower**: every hillslope solve must build a D-infinity flow graph to identify (non-existent) channels. Cidre routing adds further cost.
- **`roering_implicit` runs in this regime.** It did not before 2026-08-27: the transient hillslope step was skipped whenever the channel set was empty, on the reasoning that without Dirichlet channel sinks the solve degenerates to meaningless sink-free diffusion. That reasoning holds only for a **steady** solve with no Dirichlet data anywhere, which is singular. It does not hold for a transient step: a diffusion step over `dt` is well posed under any boundary condition and simply redistributes mass; in the coupled regime the fluvial term is itself the sink, applied to every cell in the same timestep; and an `'o'` edge supplies base level regardless of channels. Since an empty channel set is exactly what `coupled` produces on purpose, the guard was disabling the emergent-drainage-density mode it was meant to support. **The hillslope operator needs boundary conditions, not channels.** The one case still refused is `roering_implicit` + `run_steady_state` with no `'o'` edge.

### Coupling models in the literature

These two choices reflect two lineages:

| Lineage | Channel treatment | Routing | Regime name |
|---|---|---|---|
| MUDDPILE | Fixed threshold, Dirichlet | D8 | `channels_as_boundaries` |
| DonkModel | No threshold, coupled | D-infinity | `coupled` |

lsdtt3 supports both. If you want to model long-lived landscape evolution on a known network (e.g., a catchment where channels are mapped), use imposed drainage with D8. If you want to ask "what network would emerge?" or are studying **channelization**, use emergent drainage with D-infinity and Cidre.

---

## Summary: choosing a model

Here is a decision tree:

1. **Do you want to see transient dynamics** (response to climate change, base-level fall, etc.)? → **Time stepping** (`run_steady_state: false`) with one of the explicit solvers (`linear`, `roering`, or `cidre`), or `roering_implicit` for unconditional stability on fine grids.

2. **Do you want a landscape in steady state quickly?** → **Analytical steady state** (`run_steady_state: true`):
   - **Linear hillslopes:** use the SOR solver (fast, exact).
   - **Nonlinear hillslopes:** Newton is the default and needs no relaxation factor. (If you force SOR with `hillslope_force_sor`, set `hillslope_sor_omega` to 0.3–0.5 — the auto value is the linear optimum and diverges for the nonlinear operator.)
   - Optional: add snapping (`steady_state_n_cycles > 1`) or the equilibration cycle (`use_landscape_equilibration_cycle: true`) to let drainage divides migrate.

3. **Do you have a field-mapped channel network with imposed drainage?** → Set `incision_threshold_specific_area > 0`. Use D8 routing, with `roering_implicit` for transient runs.

4. **Do you want the drainage network to emerge from the landscape dynamics?** → Set `incision_threshold_specific_area: 0`, use D-infinity and Cidre (coupled routing). `roering_implicit` works here too for transient runs; give it an `'o'` edge if you also want a steady solve.

5. **How long is your run?**
   - Short transient (synthetic, < 1 My) with explicit scheme: any solver is fine.
   - Long transient (realistic timescale, > 10 My): use `roering_implicit` or `adi` for unconditional stability.
   - Steady-state: Newton is now the recommended solver for nonlinear work (fast, robust); SOR for linear work.

---

## References

Every citation below was verified against Crossref or the publisher. Where the
text names a paper, this is the paper it means.

**Fluvial incision**

- **Howard, A. D., & Kerby, G. (1983).** Channel changes in badlands.
  *Geological Society of America Bulletin*, 94(6), 739–752.
  doi:10.1130/0016-7606(1983)94&lt;739:CCIB&gt;2.0.CO;2 — the badlands study from
  which the area–slope form of the incision law comes.
- **Whipple, K. X., & Tucker, G. E. (1999).** Dynamics of the stream-power
  river incision model: implications for height limits of mountain ranges,
  landscape response timescales, and research needs. *Journal of Geophysical
  Research: Solid Earth*, 104(B8), 17661–17674. doi:10.1029/1999JB900120 —
  what the exponents `m` and `n` absorb, including runoff efficiency and
  channel-width scaling.
- **Harel, M.-A., Mudd, S. M., & Attal, M. (2016).** Global analysis of the
  stream power law parameters based on worldwide ¹⁰Be denudation rates.
  *Geomorphology*, 268, 184–196. doi:10.1016/j.geomorph.2016.05.035 — the
  empirical case that `n > 1`.
- **Braun, J., & Willett, S. D. (2013).** A very efficient O(n), implicit and
  parallel method to solve the stream power equation governing fluvial incision
  and landscape evolution. *Geomorphology*, 180–181, 170–179.
  doi:10.1016/j.geomorph.2012.10.008 — the implicit, topologically ordered
  scheme (FastScape).
- **Anand, S. K., Hooshyar, M., & Porporato, A. (2020).** Linear layout of
  multiple flow-direction networks for landscape-evolution simulations.
  *Environmental Modelling & Software*, 133, 104804.
  doi:10.1016/j.envsoft.2020.104804 — extends the O(n) ordering of Braun and
  Willett from single- to multiple-flow-direction networks, which is what makes
  a D-infinity fluvial step affordable.
- **Anand, S. K., Bonetti, S., Camporeale, C., & Porporato, A. (2022).**
  Inception of Regular Valley Spacing in Fluvial Landscapes: A Linear Stability
  Analysis. *Journal of Geophysical Research: Earth Surface*, 127(11),
  e2022JF006716. doi:10.1029/2022JF006716 — a linear stability analysis of when
  a fluvial landscape channelises, and at what spacing. Directly relevant to the
  emergent drainage density discussed above.

**Hillslope transport**

- **Culling, W. E. H. (1960).** Analytical theory of erosion. *The Journal of
  Geology*, 68(3), 336–344. doi:10.1086/626663 — introduces the analogy between
  hillslope erosion and heat conduction, i.e. linear diffusion.
- **Culling, W. E. H. (1963).** Soil creep and the development of hillside
  slopes. *The Journal of Geology*, 71(2), 127–161. doi:10.1086/626891 — derives
  the diffusive form from random particle motion. Cite alongside the 1960 paper
  when the *mechanism* rather than the equation is the point.
- **Roering, J. J., Kirchner, J. W., & Dietrich, W. E. (1999).** Evidence for
  nonlinear, diffusive sediment transport on hillslopes and implications for
  landscape morphology. *Water Resources Research*, 35(3), 853–870.
  doi:10.1029/1998WR900090 — the nonlinear flux law with a critical gradient.
- **Carretier, S., Martinod, P., Reich, M., & Goddéris, Y. (2016).** Modelling
  sediment clasts transport during landscape evolution. *Earth Surface
  Dynamics*, 4, 237–251. doi:10.5194/esurf-4-237-2016 — the Cidre model.

**Spatially varying erodibility**

- **Perne, M., Covington, M. D., Thaler, E. A., & Myre, J. M. (2017).** Steady
  state, erosional continuity, and the topography of landscapes developed in
  layered rocks. *Earth Surface Dynamics*, 5, 85–100.
  doi:10.5194/esurf-5-85-2017 — steady slopes across vertical contacts (Eq. 4)
  and horizontal layers.

**Erosivity conventions and valley spacing**

- **Perron, J. T., Dietrich, W. E., & Kirchner, J. W. (2008).** Controls on the
  spacing of first-order valleys. *Journal of Geophysical Research: Earth
  Surface*, 113, F04016. doi:10.1029/2007JF000977 — the total-area convention
  and the landscape Péclet number.
- **Perron, J. T., Kirchner, J. W., & Dietrich, W. E. (2009).** Formation of
  evenly spaced ridges and valleys. *Nature*, 460, 502–505.
  doi:10.1038/nature08174
- **Bonetti, S., Hooshyar, M., Camporeale, C., & Porporato, A. (2020).**
  Channelization cascade in landscape evolution. *Proceedings of the National
  Academy of Sciences*, 117(3), 1375–1382. doi:10.1073/pnas.1911817117 — the
  specific-area convention, adopted to remove grid dependence from the erosion
  law, and the channelisation index.
- **Litwin, D. G., et al. (2022).**
  Groundwater affects the geomorphic and hydrologic properties of coevolved
  landscapes. *Journal of Geophysical Research: Earth Surface*, 127(1),
  e2021JF006239. doi:10.1029/2021JF006239 — also uses area per unit contour
  width.
- **Theodoratos, N., Hergarten, S., & Narteau, C. (2018).** On the scaling of
  fluvial incision and hillslope diffusion in landscape evolution models.
  *Earth Surface Dynamics*, 6, 779–793. https://doi.org/10.5194/esurf-6-779-2018
  — total-area characteristic scaling and the `L_h` formulation.
- **Litwin, D. G., Malatesta, L. C., & Sklar, L. S. (2025).** Hillslope
  diffusion and channel steepness in landscape evolution models. *Earth
  Surface Dynamics*, 13(2), 277–293.
  https://doi.org/10.5194/esurf-13-277-2025 — specific-area characteristic
  scaling and the `L_g` formulation.

**Flow routing**

- **Tarboton, D. G. (1997).** A new method for the determination of flow
  directions and upslope areas in grid digital elevation models. *Water
  Resources Research*, 33(2), 309–319. doi:10.1029/96WR03137
