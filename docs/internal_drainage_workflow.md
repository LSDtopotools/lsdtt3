# Internal-drainage workflow

Internal drainage describes terrain where water terminates within the domain
rather than reaching an external boundary. Treat its routing assumptions as an
explicit modelling choice.

## Basic sequence

1. Inspect DEM validity and edge behaviour.
2. Choose boundary conditions and any intended outlet/sink locations.
3. Use depression preprocessing only when it matches the physical question.
4. Run the internal-drainage analysis and inspect its basin/depression outputs.
5. Validate results against known lakes, closed basins, or mapped drainage.

Avoid automatically filling every depression when closed basins are themselves
the feature of interest. Conversely, use a conditioned DEM when the aim is a
connected surface-flow network.

The [internal-drainage command reference](cli_reference/cli_reference_lsdtt_internal_drainage.md)
lists command-line inputs and outputs. Python users can build the equivalent
flow and basin context with `Raster`, `FlowInfo`, and the basin tools.
