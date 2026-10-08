// lsdtt_hype_routing_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// HYPE river routing. See the header. "file:line" references are to HYPE 5.36.0.

#include "lsdtt_hype_routing_cpu.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace lsdtt
{

namespace
{
constexpr double k_realzero = 1.0e-37;   // general_wc.f90:29

// remove_water (general_wc.f90:141)
void remove_water(double& vol, double q)
{
    if (vol - q > 0.0) vol = vol - q;
    else if (vol - q > -k_realzero) vol = 0.0;
    else throw std::runtime_error("hype routing: outflow exceeds the damping-box volume "
                                  "(HYPE would stop: remove water routine Q>VOL)");
}
} // namespace

hype_river make_hype_river(double length_m, double rivvel, double damp,
                           double dt_seconds, double dead_m3)
{
    hype_river r;
    r.length_m = length_m;
    r.dead_m3  = dead_m3;
    r.water_m3 = dead_m3;                                   // (:1024)
    if (rivvel == 0.0)                                      // (:1271) no regional estimation here
    {
        r.transtime = 0.0;
        r.ttstep    = 0;
        r.ttpart    = 0.0;
        r.r1 = 1.0;
        r.r2 = 1.0;
    }
    else
    {
        double vel  = rivvel;
        double part = damp;
        if (vel <= 0.001) vel = 0.001;                      // (:1284)
        if (part > 1.0) part = 1.0;
        else if (part < 0.0) part = 0.0;
        // calculate_translation_time (:803, :812-813)
        const double totaltime = length_m / vel / dt_seconds;
        r.transtime = (1.0 - part) * totaltime;
        r.ttstep    = static_cast<int>(r.transtime);          // INT truncates (:1297)
        r.ttpart    = r.transtime - static_cast<double>(r.ttstep);
        const double kt = part * totaltime;                   // (:1299)
        if (kt > 0.0)
        {
            r.r1 = 1.0 - kt + kt * std::exp(-1.0 / kt);       // (:1301)
            r.r2 = 1.0 - std::exp(-1.0 / kt);                 // (:1302)
        }
        else
        {
            r.r1 = 1.0;
            r.r2 = 1.0;
        }
    }
    r.queue.assign(static_cast<std::size_t>(r.ttstep) + 2, 0.0);
    return r;
}

double hype_river_translate(hype_river& r, double qin_m3s, double dt_seconds)
{
    // translation_in_river (sw_proc.f90:1237)
    //
    // Ring buffer: HYPE's qqueue(j) is r.queue[(r.head + j) % n], n = ttstep+2.
    // HYPE's per-step shift qqueue(1:y+1) = qqueue(0:y) (:1278) moves every
    // logical slot up by one; here it is a decrement of the head, O(1) instead
    // of O(n). The slot that wraps round to the new logical 0 holds the old
    // qqueue(y+1), where HYPE would keep the old qqueue(0); either way it is
    // overwritten by the next step's inflow before it is read, so the values
    // read (slots 1..y+1) and the arithmetic on them are unchanged.
    const std::size_t n = r.queue.size();
    r.queue[r.head] = qin_m3s * dt_seconds;                  // qqueue(0)
    if (r.transtime > 0.0)
    {
        const std::size_t y = static_cast<std::size_t>(r.ttstep);
        const double x = r.ttpart;
        std::size_t py  = r.head + y;                        // physical slot of qqueue(y)
        if (py >= n) py -= n;
        std::size_t py1 = py + 1;                            // physical slot of qqueue(y+1)
        if (py1 >= n) py1 -= n;
        double qout = (1.0 - x) * r.queue[py] + x * r.queue[py1];   // (:1268)
        qout = qout / dt_seconds;
        // queue(1:y+1) = queue(0:y) (:1278), as a head decrement.
        r.head = (r.head == 0) ? n - 1 : r.head - 1;
        return qout;
    }
    return qin_m3s;
}

double hype_river_in_transit_m3(const hype_river& r) noexcept
{
    // model_hype.f90:1424: SUM(qqueue(1:ttstep)) + qqueue(ttstep+1)*ttpart,
    // summed in HYPE's logical order (k = 1..ttstep ascending, then the
    // partial slot) through the ring-buffer mapping (see hype_river_translate).
    double v = 0.0;
    if (r.length_m > 0.0 && r.transtime > 0.0)
    {
        const std::size_t n = r.queue.size();
        const std::size_t y = static_cast<std::size_t>(r.ttstep);
        std::size_t p = r.head;                              // logical slot 0
        for (std::size_t k = 1; k <= y; ++k)
        {
            if (++p == n) p = 0;                             // logical slot k
            v += r.queue[p];
        }
        if (++p == n) p = 0;                                 // logical slot y+1
        v += r.queue[p] * r.ttpart;
    }
    return v;
}

double hype_river_step(hype_river& r, double qin_m3s, double dt_seconds)
{
    if (!(r.length_m > 0.0)) return qin_m3s;                  // no river (:6747)
    const double transq = hype_river_translate(r, qin_m3s, dt_seconds);
    r.water_m3 = r.water_m3 + transq * dt_seconds;            // (:6762; main river :2807)
    const double riverboxinflow = transq * dt_seconds;        // (:6764; main river :2809)
    // Damping box, p_rivermodel = 0 (:6881-6884; main river :2994)
    double dampq = r.r1 * riverboxinflow + r.r2 * (r.water_m3 - r.dead_m3 - riverboxinflow);
    if (dampq > r.water_m3 - r.dead_m3) dampq = r.water_m3 - r.dead_m3;
    if (dampq < 0.0) dampq = 0.0;                             // (:6887)
    remove_water(r.water_m3, dampq);
    return dampq / dt_seconds;                                // (:6894)
}

void hype_route_network_step(hype_river_network& net, const std::vector<double>& runoff_mm,
                             double dt_seconds, hype_network_flows& flows)
{
    const std::size_t n = net.local.size();
    if (runoff_mm.size() != n || net.main.size() != n || net.downstream.size() != n ||
        net.area_m2.size() != n || net.order.size() != n)
    {
        throw std::invalid_argument("hype routing: network vectors differ in length");
    }
    flows.local_outflow_m3s.assign(n, 0.0);
    flows.outflow_m3s.assign(n, 0.0);
    std::vector<double> accinflow(n, 0.0);                    // (:2610)
    for (const std::size_t s : net.order)
    {
        const double qin = runoff_mm[s] * net.area_m2[s] / (dt_seconds * 1000.0);   // (:2656)
        const double qloc = hype_river_step(net.local[s], qin, dt_seconds);
        flows.local_outflow_m3s[s] = qloc;                    // cloc (:3616)
        const double qmain_in = qloc + accinflow[s];           // (:2743)
        const double qout = hype_river_step(net.main[s], qmain_in, dt_seconds);
        flows.outflow_m3s[s] = qout;                          // cout (:3810)
        if (net.downstream[s] >= 0) accinflow[static_cast<std::size_t>(net.downstream[s])] += qout;
    }
}

} // namespace lsdtt
