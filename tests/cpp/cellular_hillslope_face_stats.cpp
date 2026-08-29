// cellular_hillslope_face_stats.cpp
//
// Research driver for the cellular hillslope plan. Reports unique face-slope
// occupancy and connected high-slope regions from an existing DEM.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_hillslope_operator_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{

struct edge
{
    std::size_t a;
    std::size_t b;
    double ratio;
};

struct disjoint_set
{
    std::vector<std::size_t> parent;
    std::vector<std::size_t> size;

    explicit disjoint_set(std::size_t n) : parent(n), size(n, 1)
    { std::iota(parent.begin(), parent.end(), 0); }

    std::size_t find(std::size_t x)
    {
        while (parent[x] != x)
        {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void join(std::size_t a, std::size_t b)
    {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (size[a] < size[b]) std::swap(a, b);
        parent[b] = a;
        size[a] += size[b];
    }
};

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 4)
    {
        std::cerr << "usage: cellular-hillslope-face-stats DEM [Sc] [NESW_BC]\n";
        return 2;
    }

    const double Sc = argc > 2 ? std::stod(argv[2]) : 0.8;
    const std::string bc = argc > 3 ? argv[3] : "opop";
    if (!(Sc > 0.0))
    {
        std::cerr << "Sc must be > 0\n";
        return 2;
    }

    lsdtt::lsdtt_data_loader loader;
    try { loader.open_raster(argv[1]); }
    catch (const std::exception& e)
    {
        std::cerr << "cannot open DEM: " << e.what() << '\n';
        return 1;
    }
    lsdtt::lsdtt_raster_cpu dem(loader);
    const auto& meta = dem.get_meta();
    const auto& mask = dem.get_mask();
    const auto dom = lsdtt::lsdtt_hillslope_domain::make(
        meta.n_rows, meta.n_cols, bc, mask);
    const lsdtt::lsdtt_hillslope_stencil<float> st(
        dom, dem.get_data().data(), meta.cell_size);

    const bool periodic_ew = dom.bc_e == 'p' && dom.bc_w == 'p';
    const bool periodic_ns = dom.bc_n == 'p' && dom.bc_s == 'p';
    std::vector<edge> edges;
    edges.reserve(dem.get_data().size() * 2);

    auto fixed = [&](int r, int c)
    { return dom.is_dirichlet_edge(r, c); };

    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
        {
            const std::size_t i = dom.index(r, c);
            if (!mask[i]) continue;

            if (c + 1 < meta.n_cols || (c == meta.n_cols - 1 && periodic_ew))
            {
                int er = r, ec = c;
                dom.neighbour(r, c, 0, 1, er, ec);
                const std::size_t ei = dom.index(er, ec);
                if (ei != i && mask[ei] && !(fixed(r, c) && fixed(er, ec)))
                {
                    const double sn = (st.at(er, ec) - st.at(r, c)) / meta.cell_size;
                    const double tang = 0.5 * (st.slope_y(r, c) + st.slope_y(er, ec));
                    edges.push_back({i, ei, std::sqrt(sn * sn + tang * tang) / Sc});
                }
            }

            if (r + 1 < meta.n_rows || (r == meta.n_rows - 1 && periodic_ns))
            {
                int sr = r, sc = c;
                dom.neighbour(r, c, 1, 0, sr, sc);
                const std::size_t si = dom.index(sr, sc);
                if (si != i && mask[si] && !(fixed(r, c) && fixed(sr, sc)))
                {
                    const double sn = (st.at(sr, sc) - st.at(r, c)) / meta.cell_size;
                    const double tang = 0.5 * (st.slope_x(r, c) + st.slope_x(sr, sc));
                    edges.push_back({i, si, std::sqrt(sn * sn + tang * tang) / Sc});
                }
            }
        }

    double max_ratio = 0.0;
    for (const edge& e : edges) max_ratio = std::max(max_ratio, e.ratio);
    std::cout << "dem," << argv[1] << '\n'
              << "rows," << meta.n_rows << '\n'
              << "cols," << meta.n_cols << '\n'
              << "active_faces," << edges.size() << '\n'
              << "max_slope_ratio," << max_ratio << '\n'
              << "threshold,n_faces,fraction,n_components,largest_component_cells\n";

    for (double threshold : {0.70, 0.80, 0.85, 1.00})
    {
        disjoint_set sets(dem.get_data().size());
        std::vector<std::uint8_t> active(dem.get_data().size(), 0u);
        std::size_t n_faces = 0;
        for (const edge& e : edges)
            if (e.ratio >= threshold)
            {
                ++n_faces;
                active[e.a] = active[e.b] = 1u;
                sets.join(e.a, e.b);
            }

        std::vector<std::uint8_t> seen(dem.get_data().size(), 0u);
        std::size_t components = 0;
        std::size_t largest = 0;
        for (std::size_t i = 0; i < active.size(); ++i)
            if (active[i])
            {
                const std::size_t root = sets.find(i);
                if (!seen[root])
                {
                    seen[root] = 1u;
                    ++components;
                    largest = std::max(largest, sets.size[root]);
                }
            }

        const double fraction = edges.empty()
            ? 0.0 : static_cast<double>(n_faces) / edges.size();
        std::cout << threshold << ',' << n_faces << ',' << fraction << ','
                  << components << ',' << largest << '\n';
    }
    return 0;
}
