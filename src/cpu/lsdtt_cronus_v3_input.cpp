// =============================================================================
// lsdtt_cronus_v3_input.cpp
//
// See lsdtt_cronus_v3_input.hpp for the format and the rationale.
// =============================================================================

#include "lsdtt_cronus_v3_input.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace lsdtt
{

namespace
{

// lsdtt3 spells nuclides "Be10"; the calculator spells them "Be-10".
std::string cronus_nuclide_name(const std::string& nuclide)
{
    if (nuclide == "Be10") return "Be-10";
    if (nuclide == "Al26") return "Al-26";
    if (nuclide == "C14")  return "C-14";
    return nuclide;   // pass through; the calculator will reject it loudly
}

// The standard the calculator assumes when the field is omitted.  C-14 has no
// standardisation field in the v3 input, so it returns empty.
std::string default_standardisation(const std::string& nuclide)
{
    if (nuclide == "Be10") return "07KNSTD";
    if (nuclide == "Al26") return "KNSTD";
    return "";
}

}  // namespace


std::string cronus_v3_sanitise_id(const std::string& raw,
                                  const std::string& fallback)
{
    std::string out;
    out.reserve(raw.size());
    for (const char c : raw)
    {
        const bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                        (c >= '0' && c <= '9') || c == '_' || c == '-';
        out += ok ? c : '_';
    }
    // A name of nothing but separators sanitises to underscores, which is legal
    // but useless as an identifier; treat it as absent.
    const bool all_underscore =
        out.empty() || out.find_first_not_of('_') == std::string::npos;
    return all_underscore ? fallback : out;
}


std::string format_cronus_v3_sample(const CronusV3Sample& sample)
{
    if (sample.nuclides.empty())
        return "";

    std::ostringstream os;
    os << std::fixed;

    const std::string id = cronus_v3_sanitise_id(sample.name, "sample");

    // Site record.  Six decimal places on the coordinates is ~0.1 m, well
    // inside any DEM's resolution; the pressure is given four because the
    // scaling responds to fractions of a hPa.
    os << id << ' '
       << std::setprecision(6) << sample.latitude_deg  << ' '
       << std::setprecision(6) << sample.longitude_deg << ' '
       << std::setprecision(4) << sample.elevation_or_pressure << ' '
       << (sample.use_pressure ? "pre" : "std") << ' '
       << std::setprecision(4) << sample.thickness_cm  << ' '
       << std::setprecision(4) << sample.density_g_cm3 << ' '
       << std::setprecision(6) << sample.shielding     << ' '
       << std::setprecision(6) << sample.erate_guess_g_cm2_yr << ' '
       << sample.sample_year << " ;\n";

    for (const CronusV3Nuclide& n : sample.nuclides)
    {
        const std::string std_name = n.standardisation.empty()
            ? default_standardisation(n.nuclide)
            : n.standardisation;

        os << id << ' ' << cronus_nuclide_name(n.nuclide) << " quartz "
           << std::setprecision(6) << n.conc_atoms_g  << ' '
           << std::setprecision(6) << (n.sigma_atoms_g > 0.0 ? n.sigma_atoms_g : 0.0);
        if (!std_name.empty())
            os << ' ' << std_name;
        os << " ;\n";
    }

    return os.str();
}


std::size_t write_cronus_v3_input_file(const std::string&                 path,
                                       const std::vector<CronusV3Sample>& samples)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open CRONUS v3 input file: " + path);

    // The calculator keys results by id, so two samples sharing one would come
    // back indistinguishable.  Suffix duplicates rather than dropping them.
    std::unordered_map<std::string, int> seen;
    std::size_t n_written = 0;

    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        CronusV3Sample s = samples[i];
        s.name = cronus_v3_sanitise_id(s.name, "S" + std::to_string(i + 1));

        const int count = ++seen[s.name];
        if (count > 1)
            s.name += "_" + std::to_string(count);

        const std::string block = format_cronus_v3_sample(s);
        if (block.empty())
            continue;
        f << block;
        ++n_written;
    }

    if (!f.good())
        throw std::runtime_error("Failed while writing CRONUS v3 input file: " + path);

    return n_written;
}

}  // namespace lsdtt
