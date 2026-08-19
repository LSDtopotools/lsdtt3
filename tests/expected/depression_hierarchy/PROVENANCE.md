# Depression-hierarchy oracle fixtures — provenance

These 25 CSV files are the **expected output of the upstream Barnes depression
hierarchy**, used by `unit_depression_hierarchy_oracle` to check lsdtt3's port in
`src/cpu/lsdtt_depression_hierarchy_cpu`.

They were generated **once, by hand, outside this repository**. lsdtt3 does not
build, link, or invoke upstream, and you do not need upstream installed to run
the test. That is deliberate: upstream pulls in the whole of RichDEM as a git
submodule plus GDAL, which is far too much to take on as a build dependency for a
test oracle.

The input DEMs are copied verbatim from upstream's `test_cases/` into
`tests/data/depression_hierarchy/` (25 ArcASCII grids, 40 KB total).

## Pinned versions

Regenerating these files requires the same two commits, or the comparison is
against a different algorithm:

| repository | commit |
|---|---|
| `github.com/r-barnes/Barnes2019-DepressionHierarchy` | `72bb3f5e92b2fd7e7f2d5627636cc407de21527b` |
| `github.com/r-barnes/richdem` (submodule `submodules/richdem`) | `af61c65584d269375a559a11ae317be3af5cd835` |

Generated 2026-08-19.

## How to regenerate

1. Clone upstream and initialise the submodule:

   ```
   git clone git@github.com:r-barnes/Barnes2019-DepressionHierarchy.git
   cd Barnes2019-DepressionHierarchy
   git checkout 72bb3f5e92b2fd7e7f2d5627636cc407de21527b
   git submodule update --init --recursive
   ```

2. Compile the dumper below against upstream's headers. It is **not** part of any
   build system; it exists only to emit the fields this test compares.
   `-DUSEGDAL` is required or RichDEM throws "RichDEM was not compiled with
   GDAL!" at run time.

   ```
   g++ -O2 -std=c++17 -fopenmp -DUSEGDAL -o dh_oracle_dump dh_oracle_dump.cpp \
       -I<upstream>/include -I<upstream>/submodules \
       -I<upstream>/submodules/richdem/include \
       $(gdal-config --cflags) $(gdal-config --libs)
   ```

3. Run it over each fixture, from the upstream root, with ocean level `-9999`:

   ```
   for f in test_cases/*.dem; do
     b=$(basename $f .dem)
     ./dh_oracle_dump $f -9999 2>/dev/null > <lsdtt3>/tests/expected/depression_hierarchy/$b.csv
   done
   ```

   RichDEM writes progress bars to stderr; only stdout is the data.

## The dumper

```cpp
#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/gdal.hpp>
#include <cstdio>
#include <string>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char** argv)
{
    if (argc != 3) { std::fprintf(stderr, "usage: %s <dem> <ocean_level>\n", argv[0]); return 2; }
    const std::string in_name = argv[1];
    const float ocean_level   = std::stod(argv[2]);

    rd::Array2D<float> topo(in_name);
    rd::Array2D<dh::dh_label_t> label(topo.width(), topo.height(), dh::NO_DEP);
    rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);

    for (unsigned int i = 0; i < label.size(); i++)
        if (topo.isNoData(i) || topo(i) == ocean_level) label(i) = dh::OCEAN;

    auto deps = dh::GetDepressionHierarchy<float, rd::Topology::D8>(topo, label, flowdirs);

    std::printf("dep,parent,lchild,rchild,pit_elev,out_elev,cell_count,ocean_parent,n_ocean_linked\n");
    for (unsigned int i = 1; i < deps.size(); i++)
    {
        const auto& d = deps[i];
        std::printf("%u,%u,%u,%u,%.9g,%.9g,%u,%d,%zu\n",
                    i, d.parent, d.lchild, d.rchild,
                    (double)d.pit_elev, (double)d.out_elev,
                    d.cell_count, d.ocean_parent ? 1 : 0,
                    d.ocean_linked.size());
    }
    return 0;
}
```

## Two things that will bite you

**`%.9g`, not `%.6g`.** Nine significant digits round-trip a `float` exactly; six
do not. An earlier version used `%.6g` and the suite failed on `testdem9` with
`1010.861` versus `1010.86` — a fixture-precision artefact that looked like a
port defect. For the same reason the test compares elevations **as `float`**:
`static_cast<double>(some_float)` does not equal `std::stod` of its `%.9g` text,
so a double comparison fails on values that are bit-identical.

**Elevations are compared exactly, with no tolerance.** They are copied from the
DEM, never computed, so any difference is a real disagreement. Do not add a
tolerance to make a failing case pass — diagnose it instead.

## Ocean semantics

Upstream's driver defines ocean as `topo.isNoData(i) || topo(i) == ocean_level`
and nothing else. lsdtt3 natively derives ocean from its boundary-condition
string, which is a different question. The test therefore sets:

```cpp
opts.ocean_from_boundary = false;   // switch off the lsdtt3 rule
opts.nodata_is_ocean     = true;
opts.ocean_level         = -9999.0f;
```

so both implementations are asked the same question with no fixture surgery.
Those options were added for exactly this reason (plan section 4.2b); they are
also useful in their own right, since a DEM whose sea is a nodata region or a
constant elevation could not previously be described at all.

## Known divergences

13 of the 25 fixtures disagree, all from **one cause**: upstream records
depressions joined through the ocean in a `Depression::ocean_linked` list and
leaves their `parent` as `OCEAN`, whereas lsdtt3 links them into the binary tree
— either chaining parents (`kerry_test3`: upstream `0,0,0`, ours `2,3,0`) or
inventing a meta-depression (`testdem2`: upstream has 2 ocean-linked leaves, ours
adds a parent node).

They are listed in `k_known_divergent` in
`tests/cpp/test_depression_hierarchy.cpp` and **reported, not skipped**. A
fixture not on that list which disagrees is a hard failure; a listed fixture that
starts agreeing is reported too, so the list gets pruned rather than rotting.

The remaining 12 agree field for field, so this is a contained defect in the
ocean-link handling, not a general disagreement about the algorithm.
