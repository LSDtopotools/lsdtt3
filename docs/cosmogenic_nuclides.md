# Cosmogenic nuclides in lsdtt3

This page consolidates the cosmogenic nuclide capability of lsdtt3: where it
came from, what physics it implements, what choices it makes that differ from
other calculators, and how those choices are tested. It is intended to be read
before using any of the three cosmogenic command-line tools, and to be citable
when reporting results produced with them.

Four points to take away if you read nothing else:

1. **lsdtt3 simulates spatially distributed cosmogenic production** and this spatial distribution
   is also used in inversion of cosmogenic concentrations into erosion rates. 
   It is the follow on from the distributed module CAIRN (Mudd et al 2016).
2. **lsdtt3 is not a CRONUS clone.** Its defaults are a physical parameter set,
   not an emulation target. Where you want agreement with the CRONUS v3 online
   calculator, there is a dedicated pair of muon schemes that deliver it to
   about a quarter of a per cent; those are opt-in.
3. **lsdtt3 applies no topographic shielding correction by default.** That is a
   deliberate methodological position, following DiBiase (2018), and it is the
   single largest reason an lsdtt3 rate will differ from an OCTOPUS rate for the
   same basin.
4. **No time-varying production correction is applied unless you ask for it.**
   Production is evaluated at the present-day field out of the box, under either
   scaling scheme.

---

## 1. Lineage: from CAIRN to lsdtt3

The method implemented here descends directly from **CAIRN** (Mudd et al., 2016,
*Earth Surface Dynamics* 4, 655–674), the basin-averaged erosion rate method
developed in LSDTopoTools. CAIRN's central idea is that a basin-averaged
cosmogenic erosion rate should be computed by evaluating production **at every
pixel** of the contributing area — with each pixel's own latitude, elevation,
pressure and shielding — and then averaging, rather than by evaluating
production once at a lumped "basin mean" position. The forward model is a
per-pixel concentration; the erosion rate is recovered by inverting the
basin-averaged concentration.

lsdtt3 inherits that structure, and inherits its numerical core almost
unchanged. The analytical solutions were checked and re-checked by hand over
several years, and the code that evaluates them is now under an explicit
numerics lock (§8.1): the values it produces cannot be altered — even by a
refactor believed to be harmless — without a test failure. This is aimed at avoiding nasty
surprises when the code is updated.

What lsdtt3 adds on top of CAIRN:

- a second scaling scheme (LSDn) evaluated per pixel, including its
  geomagnetic time integration;
- explicit atmosphere models (ERA40, NCEP, standard) rather than one hard-wired
  choice;
- time-varying production, both at a fixed sample age and solved
  self-consistently with the erosion rate;
- a transient/step-change concentration solution and a multi-nuclide path, so
  ¹⁰Be, ²⁶Al and ¹⁴C can be carried through one erosion history together;
- a 1-D column model for arbitrary erosion and burial histories;
- purpose-built CRONUS v3 emulation schemes;
- and a test suite built around external references rather than internal
  self-consistency.

What lsdtt3 deliberately does **not** inherit from CAIRN is the topographic
shielding correction (§5). This component was the most computationally expensive
component of CAIRN and DiBiase has argued it is not necessary. 

---

## 2. The production model

### 2.1 Four production pathways

Production at effective depth `z` (mass depth, g cm⁻²) is represented as a sum
of four exponentials:

```
P(z) = P0 · S_tot · Σ_i F_i · exp(−z / Γ_i)          i = 0..3
```

| i | Pathway | Γ_i (g cm⁻²) |
|---|---|---|
| 0 | neutron spallation | 160 (all schemes) |
| 1 | slow (stopped) negative muons | scheme-dependent |
| 2 | fast muons | scheme-dependent |
| 3 | "very fast" / deep muon tail | scheme-dependent |

`P0` is the sea-level high-latitude (SLHL) production rate in atoms g⁻¹ yr⁻¹,
`F_i` are the dimensionless production fractions, and `S_tot` is the total
scaling and shielding factor at the site. The `F_i` sum to 1 for most schemes;
where a published table does not sum to 1 it is **not** renormalised (the CAIRN
Table 2 ²⁶Al values sum to 0.97525 and are asserted to stay that way).

This four-exponential decomposition is the standard treatment of muon production
following Granger & Muzikar (2001) and Heisinger et al. (2002), as parameterised
by Braucher et al. (2011) and successors.

### 2.2 Steady-state concentration

For a surface eroding at a constant effective erosion rate `ε` (g cm⁻² yr⁻¹),
the concentration of a sample collected at effective depth `z` is

```
N(z) = P0 · Σ_i  F_i · Γ_i · exp(−z / Γ_i) / (ε + Γ_i · λ)
```

with `λ` the decay constant (yr⁻¹). For a sample integrated between effective
depths `z_top` and `z_bot` (a slab of finite thickness, e.g. a self-shielded
sample or a landslide scar):

```
N = P0 / (z_bot − z_top) · Σ_i F_i · Γ_i² · [exp(−z_top/Γ_i) − exp(−z_bot/Γ_i)]
                                            / (ε + Γ_i · λ)
```

Two limits are worth knowing because they are asserted in the test suite: as
`ε → 0` the concentration saturates at `N → (P0/λ)·Σ_i F_i` (finite, not a
division by zero), and for `ε ≫ Γ_0 λ` the concentration goes as `1/ε`.

### 2.3 Transient concentration

Advancing a particle by `dt` at constant `ε` from a starting effective depth `d`:

```
N(t+dt) = N(t)·exp(−λ dt)
        + P0 · Σ_i F_i · exp(−d/Γ_i) · [exp(ε dt/Γ_i) − exp(−λ dt)] / (ε/Γ_i + λ)
```

Setting `ε = 0` reduces this to pure decay plus in-growth at fixed depth,

```
N(t+dt) = N(t)·exp(−λ dt) + (P(d)/λ)·(1 − exp(−λ dt))
```

which is the standard burial expression (e.g. Lal, 1991; Dunai, 2010) and is the
form used when sediment is held at depth rather than eroded.

**Step change.** For a surface that sat at steady state under `ε_old` and then
switched to `ε_new` a time `t` ago, the sample now at the surface was at
effective depth `d = ε_new · t` at the moment of the change. So the solution is
a steady state at `ε_old` evaluated at that depth, advanced forward for `t` at
`ε_new`. No separate closed form is introduced — the two existing expressions
compose.

### 2.4 Shielding and scaling: the Vermeesch redistribution

Shielding cannot be applied as a flat multiplier on `P0`, because the pathways
attenuate at very different rates: a metre of rock removes essentially all
spallation and almost none of the deep muon flux. lsdtt3 follows Vermeesch
(2007): find the **effective depth** `x` that reproduces the required total
factor,

```
solve   Σ_i F_i · exp(−x / Γ_i)  =  S_tot        (Newton–Raphson)
then    F_i ← F_i · exp(−x / Γ_i)
```

so each pathway is attenuated by its own attenuation length rather than
uniformly. `S_tot` here is production scaling × topographic shielding, and it is
**routinely greater than 1** — Stone scaling exceeds 1 above a few hundred
metres. The solve handles that with a *negative* effective depth, amplifying
spallation most and the deep muon pathway least. (An earlier `S_tot ≥ 1` guard
silently dropped the production scaling at altitude; it was caught by the CRONUS
comparison tests, not by any internal check — see §8.3.)

Snow shielding and self shielding are handled differently: they are genuine
**mass depths** in g cm⁻², entered directly as `z_top`/`z_bot` in the equations
above, not folded into `S_tot`.

### 2.5 Muon schemes

The muon parameterisation is selected with `muon_scheme`. All five literature
schemes share the spallation attenuation length Γ₀ = 160 g cm⁻² and differ in
`P0`, in the muon attenuation lengths and in the pathway fractions.

| Scheme | P0 ¹⁰Be | P0 ²⁶Al | P0 ¹⁴C | Γ₁ | Γ₂ | Γ₃ | F(¹⁰Be) |
|---|---|---|---|---|---|---|---|
| `BraucherBorchers` **(default)** | 4.061 | 28.851 | 15.533 | 1500 | 1500 | 4320 | 0.9874 / 0.0030 / 0 / 0.0096 |
| `Braucher` | 4.30 | 31.10 | 15.21 | 1500 | 1500 | 4320 | 0.9887 / 0.0027 / 0 / 0.0086 |
| `Granger` | 5.11 | 30.31 | 5.86 | 738.6 | 2688 | 4360 | 0.9724 / 0.0186 / 0.0040 / 0.0050 |
| `Schaller` | 4.01 | 28.84 | 12.1 | 1500 | 1500 | 4320 | 0.9884 / 0.0042 / 0 / 0.0074 |
| `newCRONUS` | 4.075213 | 27.9 | 12.2 | 1459.77 | 11039.24 | 4320 | 0.98374 / 0.01372 / 0.00253 / 0 |

Decay constants: λ(¹⁰Be) = 5.00×10⁻⁷ yr⁻¹ (Chmeleff et al., 2010; Korschinek
et al., 2010), λ(²⁶Al) = 9.80×10⁻⁷ yr⁻¹. `BraucherBorchers` uses
λ(¹⁴C) = 1.2158×10⁻⁴ yr⁻¹ (5700 yr half-life); the other schemes carry the older
1.21×10⁻⁴.

Notes on provenance:

- **`BraucherBorchers`** — Braucher et al. (2011) F values with Borchers et al.
  (2016) reference production rates; ¹⁴C spallation from Lupker et al. (2015).
  This is the default because it is the most current combination of an
  independently calibrated `P0` with a physically motivated muon partition.
- **`Granger`** — the Granger & Muzikar (2001) partition, with distinct slow and
  fast muon attenuation lengths rather than the merged 1500 g cm⁻² of the
  Braucher family.
- **`newCRONUS`** — the updated ¹⁰Be production curve from the CRONUScalc line
  (Marrero et al., 2016), with the other nuclides left at COSMOCALC values.
- **`Schaller`** — Schaller et al. (2016), differing from `Braucher` mainly in
  the negative-muon fraction.

### 2.6 The CRONUS v3 emulator schemes

Two further schemes exist for one purpose only: **reproducing the CRONUS v3
online calculator's erosion rates.**

```
muon_scheme: cronusEmulatorSt      # pair with scaling_scheme: stone
muon_scheme: cronusEmulatorLSDn    # pair with scaling_scheme: lsdn
```

CRONUS v3 ships two muon models — a single-exponential "fast approximation" with
a pressure-dependent effective attenuation length, and a numerically integrated
model-1A muon inventory tabulated on a (pressure, log erosion rate) grid. **Its
erosion-rate root-find uses the table; the exponential is used only for exposure
ages, saturation checks and plotting.** That distinction matters: the
exponential's attenuation length is a *near-surface* value, whereas an eroding
column integrates muon production to a depth of order `ε/λ` — kilometres of
rock, dominated by the deep fast-muon tail. Measured against v3's own table, the
exponential delivers only 17.6–78.5 % of the true muon inventory across the
erosion rates relevant to erosion-rate work.

lsdtt3 therefore fits the **tabulated** inventory with the three muon pathways it
already has:

```
N_mu(ε, p) = Σ_k  P_k(p) · Λ_k(p) / (ε + Λ_k(p) · λ)          k = 1..3
```

with `log P_k` and `log Λ_k` quadratic in pressure `p`, and `Λ_1 < Λ_2 < Λ_3`
enforced structurally so the ordering cannot swap with pressure. Residuals
against the table over 450–1020 hPa and ε from 2.7×10⁻⁶ to 2.7 g cm⁻² yr⁻¹:

| exponentials used | ¹⁰Be | ²⁶Al | ¹⁴C |
|---|---|---|---|
| 1 (a single-exponential scheme) | 22.4 % | 19.8 % | 10.6 % |
| 2 | 2.8 % | 2.9 % | 1.9 % |
| **3 (shipped)** | **0.77 %** | **0.87 %** | **1.15 %** |

The surface muon production against pressure follows Balco (2017, Eq. 7),

```
P_mu,0(p) = P_mu,SLHL · exp((1013.25 − p) / L_mu,atm)
```

verified against the v3 table in the `ε → 0` limit to 1–3 % over 700–1013 hPa —
though the fit is deliberately **not** anchored to that formula, because the
table itself departs from it by up to 3 % and anchoring measurably degrades the
fit.

Three practical constraints:

1. These schemes are **pressure- and nuclide-dependent**, so a pressure raster
   is computed automatically and one parameter set serves one nuclide.
2. They are **steady-state only**; the step-change and transient paths reject
   them explicitly.
3. `lsdtt-cosmo-column-tool` does not offer them — it works on a single column
   with no pressure raster.

---

## 3. Scaling schemes

`scaling_scheme` selects how spallation production is scaled from SLHL to the
site. The muon pathways are unaffected by this choice directly; they receive the
site factor through the Vermeesch redistribution of §2.4.

### 3.1 `stone` — Stone (2000) / Lal (1991)

The default, and the scheme CAIRN and OCTOPUS use. A latitude- and
pressure-dependent scaling factor built from Stone's tabulated coefficients,
with the spallation fraction `Fsp` separating the neutron and muon components of
the original formulation. Implemented as a faithful port of the LSDTopoTools2
routine (of which CAIRN is a part), itself following Balco et al. (2008).

### 3.2 `lsdn` — Lifton, Sato & Dunai (2014)

Nuclide-specific, spectrum-based scaling, interpolated from the LSDn2015 grids
shipped in `crn_parameters.nc`. Supports ¹⁰Be, ²⁶Al and ¹⁴C only. The grids give
coefficients `b(Rc, p)` and `m(Rc, p)` on a cutoff-rigidity × pressure lattice,
combined with a solar modulation parameter:

```
A(Rc, p, S) = b(Rc, p) + S · m(Rc, p)
```

normalised to its SLHL value. lsdtt3 evaluates this **per pixel**, with each
pixel's own latitude, longitude and pressure.

The longitude dependence is not incidental. Under a time-varying field the
geomagnetic latitude is the great-circle distance to a *moving* palaeopole, so
the cutoff rigidity at a site depends on its longitude as well as its latitude.
A latitude-only implementation would carry an error growing with distance from
the palaeopole meridian — and a latitude-only sensitivity test would never
detect it.

### 3.3 Atmosphere

Pressure is the dominant control on production, so the atmosphere model is a
first-order choice, not a detail.

| `atmosphere` | Description |
|---|---|
| `ERA40` **(default)** | Balco/Lifton ERA40 climatology with a latitude-dependent lapse rate (COSPAR CIRA-86 polynomial) |
| `NCEP` | Balco's `NCEPatm_2` reanalysis climatology — **the atmosphere CAIRN and OCTOPUS use** |
| `standard` | SLHL barometric formula, no climatology |

The standard atmosphere is

```
P(h) = 1013.25 · exp[ (g M / R dT/dz) · ( ln T_slp − ln (T_slp − dT/dz·h) ) ]
```

with `T_slp = 288.15 K` and `dT/dz = 0.0065 K m⁻¹`. `ERA40` and `NCEP` fall back
to `standard` if their grids are absent from `crn_parameters.nc`.

**Match the atmosphere between tools.** The atmosphere sets the pressure at every
pixel, so forward-modelling under one and inverting under another biases the
recovered erosion rate silently. The ERA40-versus-NCEP difference is worth about
−0.64 % at basin scale in the OCTOPUS comparison (§6).

---

## 4. Time-varying production

**Nothing time-varying happens by default.** Out of the box, production is
evaluated at the present-day field under either scaling scheme, and every run
reports which it did:

```
Palaeomag     : off (modern field)
```

There are two ways to switch it on, and they mean different things.

### 4.1 Fixed sample age

`sample_age_yr` sets an exposure age in years BP and integrates production over
that window. These cosmogenic routines are for eroding landscapes, see section
on "self-consistent age" below what age means in this context (basically it is a residence time). 
What `sample_age_yr` selects depends on the scaling scheme:

- under `stone`, the Muscheler VDM palaeomagnetic correction — a time average of
  the Stone factor over the virtual dipole moment record, latitude-dependent;
- under `lsdn`, a duration-weighted average of LSDn production over the
  Pavón-Carrasco geomagnetic record — a Holocene moving palaeopole followed by a
  GAD dipole record extending to ~20 Ma. This is **longitude-dependent as well
  as latitude-dependent**, for the reason given in §3.2.

"Now" is taken as 2010, the reference year of the packed record. lsdtt3 has no
sample-collection-year input, so no (collection year − 2010) shift is applied;
this is negligible against 10²–10⁷ yr records but is an assumption rather than an
exact treatment.

### 4.2 Self-consistent age

For an eroding surface there is no independent exposure age — the integration
timescale is set by the erosion rate itself, `t = Γ_sp/ε`. `self_consistent_age:
true` solves the two together: invert for `ε`, form `t = Γ_sp/ε`, rebuild the
production raster, re-invert, iterate to a relative tolerance on `ε`.

This is the CRONUScalc / Balco v3 convention, and **it is how most published
LSDn rates were computed**. It is off by default because it changes what the
reported rate means and because it costs a production-raster rebuild per
iteration per sample. A sample that fails to converge within the iteration cap
is reported with a warning and flagged in the output CSV; its last iterate is
written, never silently accepted.

Leaving it off while using `cronusEmulatorLSDn` raises an advisory warning,
because that combination cannot reproduce CRONUS v3 — measured bias at one
fixture is +10.9 % with the iteration off against +0.5 % with it on.

### 4.3 How it is made fast

Evaluated naively, the time integration costs up to ~1800 grid evaluations per
pixel (Pavón-Carrasco has 91 Holocene steps plus up to 1710 long-record steps).
Both schemes therefore precompute a time-averaged surface once per run and
evaluate it per pixel by interpolation.

The important property is that **the approximation is horizontal only**. For
LSDn, `b` and `m` are linear in pressure between adjacent grid rows and time
averaging is a linear operator, so building on the 380 native pressure nodes and
interpolating linearly in pressure afterwards is identical to time-averaging at
arbitrary pressure. For the Stone palaeomagnetic correction, the factor is an
interpolation over a table whose breakpoints are pressure-independent, so the
*weights* can be time-averaged and the pressure applied afterwards at full
accuracy.

The shipped horizontal resolutions were chosen by measurement, not assumed:
0.01° for the Stone surface (max relative deviation 7.9×10⁻⁸ from the direct
definition, at the precision of the float32 production raster it feeds) and 0.1°
for LSDn, whose build carries the extra factor of 380 pressure nodes.

---

## 5. Topographic shielding: why it defaults to off

lsdtt3 applies **no topographic shielding correction by default**
(`topo_shield_value: 1`). This is a deliberate departure from CAIRN and from
most older basin-averaged workflows, and it is the single largest methodological
difference from OCTOPUS.

The reasoning follows **DiBiase (2018)**. The conventional correction accounts
only for the reduction of production *at the surface*, while assuming the
effective mass attenuation length is spatially uniform. It is not: attenuation
length increases on steep slopes, and the two effects largely offset. Correcting
for one and not the other introduces a **slope-dependent systematic error** —
which is worse than a simple bias, because slope correlates with erosion rate,
so it can manufacture spurious topography–erosion relationships. In steep
landscapes the surface-only correction reaches ~20 %.

DiBiase's conclusion is that topographic shielding corrections are inappropriate
for catchment-mean erosion rates in most cases, and are needed only for steep
catchments with a **non-uniform distribution of quartz and/or erosion rate**.

Set `topo_shield_value` yourself when you are in that case, or when reproducing
a rate computed the older way. Two things to know if you do: it is applied as a
single uniform value across the whole basin, not per pixel; and in lsdtt3's
published-dataset comparisons the topographic shielding term decomposes at about
**+0.70 %** at basin scale.

The same argument does not apply to `self_shield_eff_thickness` or
`snow_shield_eff_thickness`. Those are ordinary depth-attenuation terms in
g cm⁻² and default to 0 simply because the common case is an unshielded surface
sample.

---

## 6. Why lsdtt3 and OCTOPUS give different numbers

**OCTOPUS** (Codilean et al., 2018) is a global database of basin-averaged
cosmogenic nuclide data in which erosion rates are recomputed consistently using
CAIRN. Because lsdtt3 descends from CAIRN, it is a natural comparison.

### 6.1 OCTOPUS's conventions, and lsdtt3's defaults against them

| Convention | OCTOPUS / CAIRN | lsdtt3 default | Effect |
|---|---|---|---|
| Scaling | Stone (2000)/Lal | `stone` — same | — |
| Atmosphere | NCEP (`NCEPatm_2`) | `ERA40` | ≈ −0.64 % |
| Muon scheme | CAIRN's `Braucher` | `BraucherBorchers` | scheme-dependent |
| Palaeomagnetic correction | none | none (`sample_age_yr: 0`) | — |
| Topographic shielding | applied | **not applied** | ≈ +0.70 % |
| Production evaluation | per pixel (CAIRN) | per pixel | — |
| DEM | SRTM 90 m throughout | whatever you supply | resolution term |

To run in the CAIRN convention, set all of these explicitly:

```
muon_scheme:         Braucher
scaling_scheme:      stone
atmosphere:          NCEP
sample_age_yr:       0
self_consistent_age: false
topo_shield_value:   <your basin shielding factor>
```

Note that `muon_scheme` and `atmosphere` both differ from the lsdtt3 defaults
and both change the answer on their own.

### 6.2 The measured difference, decomposed

Against three OCTOPUS-derived datasets, compared against their **LSDn** erosion
rates and run age-matched (`self_consistent_age` on, so both sides integrate
over the same erosional timescale), with the sign
`100·(lsdtt3 − published)/published`:

| Dataset | n | median | mean | min | max |
|---|---|---|---|---|---|
| Clapp 2002 (Yuma Wash) | 21 | −4.59 | −4.78 | −8.27 | −3.57 |
| Dethier 2014 (Front Range) | 12 | −5.85 | −5.85 | −8.08 | −3.71 |
| DiBiase 2010 (San Gabriel) | 14 | −8.23 | −8.76 | −21.22 | −2.65 |

**These are bands to watch, not accuracy figures.** Each mixes the production
kernel difference, the per-pixel-versus-lumped term, the atmosphere difference,
the shielding convention and DEM/delineation differences.

The decomposition was *predicted in advance* and then tested: kernel −6.14,
per-pixel versus lumped −0.53, ERA40-rather-than-NCEP −0.64, topographic
shielding +0.70 ⇒ ≈ **−6.6 %**. The per-pixel-versus-lumped term belongs in this
particular decomposition because the comparison is against the LSDn rates, which
are the lumped ones (§6.3); it would **not** appear in a comparison against
OCTOPUS's CAIRN-derived rates, which are per pixel. Re-running age-matched landed on that for two of
the three datasets. DiBiase sits wider, as it does in every comparison, because
those basins erode up to 10× faster and the spread tracks erosion rate — the
signature of a muon-scheme difference rather than a positional one.

The tier originally ran at the default `sample_age_yr: 0` while the published
values are time-averaged over Γ/ε, and measured **+3.8 to +10.1 %** — positive.
Reading that as "accuracy" would have been wrong in both magnitude and sign.
This is the clearest illustration of why the time-integration convention has to
be matched before any comparison means anything.

### 6.3 Per-pixel versus lumped: a result, not an error

**OCTOPUS's headline erosion rates are computed per pixel, using CAIRN**, in the
same spatially distributed way lsdtt3 computes them. That is not where the two
differ, and it is why the comparison in §6.2 is a meaningful one at all.

Where a lumped evaluation does appear is in OCTOPUS's **LSDn-equivalent rate**,
supplied for comparison: that one is computed from a single lumped basin
parameter — a basin-representative latitude and pressure — rather than pixel by
pixel. RIVERSAND likewise evaluates the scaling once at a lumped basin position.

So the per-pixel-versus-lumped difference is a real quantity, but it is scoped:
it applies to the **LSDn comparison rates**, not to the CAIRN-derived rates that
OCTOPUS publishes as its own. Measured over 47 basins at 5 kyr under LSDn,
comparing lsdtt3's per-pixel LSDn mean against LSDn evaluated once at the
published basin centre: median **−0.525 %**, IQR [−1.590, +0.141]. It is a
scientific result, not an error term.

That sits on top of an independent published estimate of the same quantity —
Codilean & Munack (2025, *Geochronology* 7:113) report an LSDn median of −1.0 %
with IQR −1.8 to −0.1 % over 290 points: same sign, same order of magnitude,
overlapping IQRs.

A second result falls out of the same machinery: inverting lsdtt3's per-pixel
**Stone** basin mean for the pressure that would reproduce it at the published
centre recovers OCTOPUS's own published pressure to a median of **+0.73 hPa**
(IQR −0.18 to +2.19), about 0.1 % of a ~900 hPa pressure. That confirms our
reading of their reduction.

---

## 7. Why lsdtt3 and CRONUS give different numbers

The answers are **expected to differ**, and each comparison is recorded as a
percent difference held inside a frozen band rather than asserted to be zero.

| Comparison | Difference | Cause |
|---|---|---|
| vs CRONUS-2.2 (St) | up to ~40 % in concentration, a function of `S_tot` | higher P0 (4.49) plus Heisinger muon overproduction |
| vs CRONUScalc | 2.5–5 % in ε | muon scheme plus spatial averaging |
| vs CRONUScalc running `newCRONUS` | < 2 % | residual spatial averaging and fixed Λ₀ = 160 |
| vs CRONUS v3 St, default `BraucherBorchers` | +1.8 % to +6.65 %, mean +3.6 % | Balco (2017) model-1A muons; grows with ε as the muon share rises |
| **vs CRONUS v3 St, `cronusEmulatorSt`** | **mean −0.23 %, max 1.03 % (¹⁰Be)** | purpose-built to reproduce v3 |
| vs CRONUS v3 LSDn | mean +11.8 % against time-independent rates | LSD framework and time integration |

Validation of the emulator against the **live calculator**: 210 points, latitude
{−60, −30, −1, 30, 60} × elevation {0, 500, 1000, 1500, 2500, 3500} m × ε
{0.001 … 1.0} g cm⁻² yr⁻¹, sign `100·(lsdtt3 − v3)/v3`:

| Run | mean | max abs | latitude trend |
|---|---|---|---|
| `cronusEmulatorSt`, ¹⁰Be | **−0.233 %** | 1.028 % | 0.114 pp |
| `cronusEmulatorSt`, ²⁶Al | −0.255 % | 1.083 % | 0.108 pp |
| `cronusEmulatorSt`, ¹⁴C | −0.252 % | 4.234 % | 0.255 pp |
| `BraucherBorchers`, ¹⁰Be (baseline) | −4.638 % | 10.534 % | 2.767 pp |
| `cronusEmulatorLSDn`, ¹⁰Be | −3.366 % | 26.573 % | 9.154 pp |

¹⁴C's 4.2 % worst case is confined to ε = 0.001, where the nuclide is saturated
and the inversion is ill-conditioned: a sub-percent concentration difference maps
to a several-percent erosion rate. At every other rate ¹⁴C stays inside 1.25 %.

**`cronusEmulatorLSDn` does not reach the same standard, and the cause is not
the muons.** v3's LSDn root-find integrates a time-varying spallation production
history; the residual is −8 % at both poles, +1 % at the equator and asymmetric
between −30° and +30°, which is spallation- and field-side structure. Use
`self_consistent_age: true` with it — that is what closes most of the gap.

### 7.1 A note on riversand

**riversand** is an open-source Python package that computes basin-averaged
erosion rates by submitting to the CRONUS v3 online calculator. It appears here
twice.

First, as a **route to the reference values**: CRONUS v3's MATLAB source is not
public, but the calculator is scriptable, and lsdtt3's reference CSV is fetched
by the same HTTP mechanism riversand uses. That fetch is a manual, reviewable
event — the frozen CSV is what CI reads, and **CI never touches the network**.

Second, as an **independent cross-check tool**: a harness runs riversand over the
same published datasets lsdtt3 is compared against and diffs the two summaries.
It is not part of CI (it needs downloads and a network) and it is run manually.

For interpreting the numbers, the important thing about riversand is that it
evaluates the scaling **once at a lumped basin position**. lsdtt3 evaluates it
per pixel, so a riversand-versus-lsdtt3 difference contains the
per-pixel-versus-lumped term of §6.3 by construction and should not be read as a
kernel discrepancy. Note this is *not* a property riversand shares with
OCTOPUS's published rates, which are computed per pixel by CAIRN — see §6.3.

---

## 8. Testing

The cosmogenic code is the most heavily tested part of lsdtt3. Full detail is in
[Benchmarking and testing](testing_and_benchmarking.md); this is the summary.

### 8.1 The numerics lock

The analytical solutions underpinning CAIRN and the muon-based step change were
checked by hand over several years, and the suite is built to make altering them
unnoticed impossible.

- **Sentinel pinning.** A matrix of exact concentration and erosion-rate values —
  every nuclide × scheme at several points in (ε, depth, shielding, time) — frozen
  to full double precision and compared at 1×10⁻¹² relative tolerance. This is far
  stricter than physical accuracy requires; its job is to trip on *any* change,
  including a refactor or an algebraic rearrangement. LSDn has its own separate
  sentinel so a change to one data path cannot be masked by the other. Both are
  pinned by SHA-256.
- **Exact parameter pins.** `F`, `Γ`, `P0` and `λ` for every muon scheme are
  asserted with exact equality on the doubles, not a tolerance.
- **Gated regeneration.** There is no regenerate-on-failure convenience path. A
  changed value must be reviewed and approved before it is committed.

### 8.2 Equation fidelity against the papers

The tests do not check the code against itself. The test driver carries its own
independent double-precision re-implementations of the published equations,
written from the papers, so agreement means something. Asserted properties
include: the depth-integrated form reducing to the point form as the slab
shrinks; step-change end members (`t = 0` → steady state at `ε_old`, `t → ∞` →
steady state at `ε_new`, monotone approach, trajectory bounded by the two);
Vermeesch's sum identity `Σ F_i S_i = S_tot` with the correct ordering in both
the attenuating and amplifying cases; saturation and `1/ε` limits; forward →
inverse round-trip recovery of ε across all nuclides and schemes; flux-weighted
versus unweighted basin averaging identities; and AMS standardisation factors
against the CRONUS table.

> **Internal cross-checks are not independent oracles.** The analytic basin
> formula and the particle-column integrator agree to 9.25×10⁻⁸ — but both share
> the same scaling code, so both were blind to the `S_tot ≥ 1` bug. The external
> CRONUS comparison is what caught it.

### 8.3 What the external comparisons caught

Two defects were found by comparison rather than by unit tests, which is the
argument for keeping external references in the suite:

1. **The `S_tot ≥ 1` guard** that silently discarded production scaling above a
   few hundred metres. Invisible to every internal check; exposed by the
   elevation axis of the CRONUS comparison fixture.
2. **A standardisation double-count** in the OCTOPUS comparison tier.

### 8.4 Fixtures and tiers

| Fixture | What it is |
|---|---|
| Hypothetical point site | 70° N, 1007 hPa, no shielding, ε sweep — the CAIRN paper's own comparison fixture, so the forward curve is directly comparable to its published difference curves |
| Synthetic basin | A small V-valley DEM with a single forced outlet, generated in C++, with uniform-ε and two-rate variants — exercises basin averaging and flux weighting |
| Real-coordinate sample set | 15 synthetic samples over 5 latitudes × 3 elevations × an ε cycle, both nuclides — maps the difference surface against CRONUS |
| Real-landscape tier | Published datasets downloaded, run, and deleted; compared against frozen summaries |

Comparisons publish percent differences in **more than one sign convention**
(kernel-oracle, per-pixel-versus-lumped and real-landscape tiers each run a
different way round). Figures from one convention must never be added to figures
from another; read the emitting test before combining any two numbers.

---

## 9. The three command-line tools

### 9.1 `lsdtt-cosmo-tool` — measured concentrations → erosion rates

The inverse tool, and the direct descendant of CAIRN. Loads a DEM, builds D8
flow routing, computes a production raster, snaps sample points to the channel
network, and inverts each basin-averaged concentration to an effective erosion
rate with uncertainty bounds.

Pipeline: load DEM (projected CRS required) → pit fill → flow routing →
pressure raster → production raster (Stone or LSDn) → load and snap sample
points → per sample: apply AMS standardisation, invert, propagate uncertainty →
write `<write_fname>_CN_ERATE.csv`.

Use it when you have measured nuclide concentrations and want erosion rates.

[Full parameter reference](cli_reference/cli_reference_lsdtt_cosmo_tool.md)

### 9.2 `lsdtt-cosmo-prediction` — erosion scenarios → predicted concentrations

The forward tool. Predicts concentrations across a DEM under an imposed erosion
scenario, in five modes:

| Mode | What it does |
|---|---|
| `calculate_crn_concentration_raster` | Steady-state concentration at every pixel from an erosion-rate field |
| `calculate_step_change_crn_concentration_raster` | Transient concentration after a step change in erosion rate, with a per-pixel time-since-change |
| `calculate_multi_isotope_step_change` | ¹⁰Be, ²⁶Al and ¹⁴C through the **same** step change in one particle evolution, plus the two ratio rasters |
| `calculate_multi_isotope_transient` | The same three nuclides through an arbitrary erosion history from a CSV |
| `calculate_accumulated_crn_concentration` | Erosion-weighted downstream accumulation, `Σ(C·ε)/Σ(ε)` — what a sediment sample would contain |

The multi-isotope modes exist because the **ratios** carry the transient signal:
half-lives differ by orders of magnitude (¹⁴C 5.7 kyr, ²⁶Al 0.7 Myr, ¹⁰Be
1.39 Myr), so a nuclide pair dates a perturbation in a way one nuclide cannot. A
ratio is only meaningful if both members saw an identical scenario, so all three
are computed from one particle evolution per pixel — same shielding, same rates,
same elapsed time — rather than from three separate runs. Each nuclide is
bit-for-bit what the single-nuclide mode returns.

Sampling the accumulated raster at lat/lon points writes a CSV that
`lsdtt-cosmo-tool` reads directly, closing a forward → inverse round trip
against the imposed erosion rate.

Use it to design sampling campaigns, to test what a scenario would look like in
the data, and to check inversion machinery against a known answer.

[Full parameter reference](cli_reference/cli_reference_lsdtt_cosmo_prediction.md)

### 9.3 `lsdtt-cosmo-column-tool` — one column, arbitrary history

A purely 1-D particle calculation: no DEM, no raster. Seven modes:
`steady_state`, `step_change`, `full_transient` (a column of particles driven by
an erosion-history CSV, writing a depth profile), plus `burial_pair`,
`step_pair` and their batch variants for two-nuclide inversions.

Use it for depth profiles, burial dating, sensitivity work, and for
understanding what a given erosion or burial history does to a nuclide pair
before committing to a landscape-scale run.

One convention to get right: in an erosion-history CSV, **row 0 is the most
recent** and years BP increase with index. A reversed record is rejected
fatally — before that check existed it returned a confidently wrong number with
no warning.

[Full parameter reference](cli_reference/cli_reference_lsdtt_cosmo_column_tool.md)

---

## 10. Citing and reporting

When reporting lsdtt3 cosmogenic results, state at minimum: the nuclide, the
muon scheme, the scaling scheme, the atmosphere, whether a time-integration was
applied (fixed age or self-consistent), the topographic shielding value, and the
DEM and its resolution. Any of these changes the number, and several of them
differ from the defaults of other calculators.

---

## 11. References

**Method and implementation**

- Mudd, S. M., Harel, M.-A., Hurst, M. D., Grieve, S. W. D. & Marrero, S. M.
  (2016). The CAIRN method: automated, reproducible calculation of catchment-averaged
  denudation rates from cosmogenic nuclide concentrations. *Earth Surface
  Dynamics*, 4, 655–674.
- Balco, G., Stone, J. O., Lifton, N. A. & Dunai, T. J. (2008). A complete and
  easily accessible means of calculating surface exposure ages or erosion rates
  from ¹⁰Be and ²⁶Al measurements. *Quaternary Geochronology*, 3, 174–195.

**Production, scaling and atmosphere**

- Lal, D. (1991). Cosmic ray labeling of erosion surfaces: in situ nuclide
  production rates and erosion models. *Earth and Planetary Science Letters*,
  104, 424–439.
- Stone, J. O. (2000). Air pressure and cosmogenic isotope production. *Journal
  of Geophysical Research*, 105(B10), 23753–23759.
- Lifton, N., Sato, T. & Dunai, T. J. (2014). Scaling in situ cosmogenic nuclide
  production rates using analytical approximations to atmospheric cosmic-ray
  fluxes. *Earth and Planetary Science Letters*, 386, 149–160.
- Borchers, B. et al. (2016). Geological calibration of spallation production
  rates in the CRONUS-Earth project. *Quaternary Geochronology*, 31, 188–198.

**Muons**

- Heisinger, B. et al. (2002). Production of selected cosmogenic radionuclides
  by muons: 1. Fast muons; 2. Capture of negative muons. *Earth and Planetary
  Science Letters*, 200, 345–355 and 357–369.
- Granger, D. E. & Muzikar, P. F. (2001). Dating sediment burial with in situ-produced
  cosmogenic nuclides: theory, techniques, and limitations. *Earth and Planetary
  Science Letters*, 188, 269–281.
- Braucher, R., Merchel, S., Borgomano, J. & Bourlès, D. L. (2011). Production of
  cosmogenic radionuclides at great depth: a multi element approach. *Earth and
  Planetary Science Letters*, 309, 1–9.
- Balco, G. (2017). Production rate calculations for cosmic-ray-muon-produced
  ¹⁰Be and ²⁶Al benchmarked against geological calibration data. *Quaternary
  Geochronology*, 39, 150–173.
- Schaller, M. et al. (2009, 2016). Source of the revised negative-muon
  attenuation length used by the `Schaller` scheme.
- Lupker, M. et al. (2015). *Quaternary Geochronology*. Source of the in-situ
  ¹⁴C spallation production rate used by `BraucherBorchers`.

**Shielding**

- Vermeesch, P. (2007). CosmoCalc: an Excel add-in for cosmogenic nuclide
  calculations. *Geochemistry, Geophysics, Geosystems*, 8, Q08003.
- DiBiase, R. A. (2018). Short communication: Increasing vertical attenuation
  length of cosmogenic nuclide production on steep slopes negates topographic
  shielding corrections for catchment erosion rates. *Earth Surface Dynamics*,
  6, 923–931.

**Decay constants and standards**

- Chmeleff, J., von Blanckenburg, F., Kossert, K. & Jakob, D. (2010). Determination
  of the ¹⁰Be half-life by multicollector ICP-MS and liquid scintillation counting.
  *Nuclear Instruments and Methods in Physics Research B*, 268, 192–199.
- Korschinek, G. et al. (2010). A new value for the half-life of ¹⁰Be by
  heavy-ion elastic recoil detection and liquid scintillation counting. *Nuclear
  Instruments and Methods in Physics Research B*, 268, 187–191.
- Nishiizumi, K. et al. (2007). Absolute calibration of ¹⁰Be AMS standards.
  *Nuclear Instruments and Methods in Physics Research B*, 258, 403–413.

**Calculators and databases**

- Marrero, S. M. et al. (2016). Cosmogenic nuclide systematics and the
  CRONUScalc program. *Quaternary Geochronology*, 31, 160–187.
- Codilean, A. T. et al. (2018). OCTOPUS: an open cosmogenic isotope and
  luminescence database. *Earth System Science Data*, 10, 2123–2139.
- Codilean, A. T. & Munack, H. (2025). *Geochronology*, 7, 113.
