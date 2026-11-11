# Design and algorithms

The pages in this section explain what parts of lsdtt3 do and how to use them.
Each one takes a topic, such as flow routing, chi analysis or downloading
data, and covers:

- **the method:** the equations, the algorithms, the choices made in the
  code, and the papers they come from;
- **how to use it:** which tools and parameters run it, in what order, and
  what the outputs mean.

Read them when you need to know why a tool gives the result it does, how to
choose its settings, or what to cite.

**This section is not comprehensive.** It covers the parts of lsdtt3 where the
method needs explaining; many tools and options have no page here. For what
every tool can do, use the
[command-line reference](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md):
one page per tool, listing every parameter with its default and meaning. Every
tool also prints the same information itself: `<tool> -h` lists its
parameters, and `<tool> -h <parameter>` explains one. The design and
algorithms pages are for the details behind them.

For the same analyses from Python, see the [Python guides](../python/index.md).
The Python functions run the same C++ routines as the command-line tools, so
the methods described here apply to both.
