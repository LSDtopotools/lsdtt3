# LSDTopoTools3

<p align="center">
  <img src="LSD-logo.png" alt="LSDTopoTools3 logo" width="220">
</p>

**LSDTopoTools3** (lsdtt3) is a C++17 library and set of command-line tools for topographic analysis of digital elevation models (DEMs). It is a ground-up redesign of [LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2), retaining the core algorithms while rebuilding the architecture for modern C++, optional GPU support, and easier incremental development.

---

## Quick links

| I want to… | Go to |
|---|---|
| Understand what lsdtt3 does | [Overview](overview.md) |
| Install lsdtt3 on Windows | [Windows (WSL2)](developer_notes/lsdtt3_in_windows.md) |
| Install lsdtt3 on macOS | [macOS (Multipass)](developer_notes/lsdtt3_in_macos.md) |
| Install lsdtt3 on Linux | [Build from source](installing_locally.md#build-from-source-on-linux) |
| Install pre-built executables (beta) | [Beta testing guide](developer_notes/beta_testing_notes.md) |
| Look up parameters for a tool | [Command-line reference](cli_reference/cli_reference_lsdtt_channel_extraction.md) |
| Contribute or modify the code | [Developer notes](developer_notes/lsdtt3_architecture.md) |

---

## Capabilities at a glance

**DEM preparation** — reproject, fill, breach depressions, clip to valid extent.

**Surface analysis** — hillshade, curvature, roughness, local relief, Frangi filter, cloth-drape incision depth.

**Flow routing** — D8, D-infinity, and MFD routing; topographic wetness index; 64-bit pixel indexing for very large DEMs.

**Channel networks** — channel extraction, junction networks, Strahler ordering; export to CSV, GeoPackage, or FlatGeobuf.

**Basin analysis** — basin delineation, hydromorphic statistics, chi analysis and k_sn inversion.

**Cosmogenic nuclides** — basin-averaged CRN erosion rates, concentration prediction, 1-D transient column simulation.

**Landscape evolution** — implicit stream-power incision, linear and nonlinear hillslope transport, uniform uplift.

---

## Citing lsdtt3

If you use lsdtt3 in research, please also cite the original LSDTopoTools2:

> Mudd, S. M., Clubb, F. J., Grieve, S. W. D., Milodowski, D. T., Gailleton, B., Hurst, M. D., Valters, D. V., Wickert, A. D., & Hutton, E. W. H. (2023). LSDtopotools/LSDTopoTools2: LSDTopoTools2 v0.9. Zenodo. [https://doi.org/10.5281/zenodo.8076231](https://doi.org/10.5281/zenodo.8076231)

**Copyright:** Simon M. Mudd | **Licence:** GPL-3.0
